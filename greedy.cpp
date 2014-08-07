// AUTHOR :  Scott Field 
//           sfield@umd.edu
//
// DATE: Jan 20, 2013
//
// PURPOSE: mpi version of greedy (pivoted MGS) algorithm 


//#include "nr3.h"
#include "gauss_wgts.h"
#include <libconfig.h++>
#include <mpi.h>

#include "hdf5.h"
#define FILE_H5 "file.h5"

//#include <boost/numeric/ublas/vector.hpp>
//#include <boost/numeric/ublas/io.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <complex>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_complex.h>
#include <gsl/gsl_block_complex_float.h>
//#include <gsl/gsl_errno.h>

#include "spa_waveforms.h"
#include "training_set.hpp"
#include "gsl_helper_functions.h"
#include "quadratures.h"


// *** ONLY MODEL SPECIFIC PART OF THE CODE *** //
void FillTrainingSet(gsl_matrix_complex *TS_gsl, const gsl_vector *xQuad, const gsl_vector_complex *wQuad, const TrainSet ts, const int rank)
{

    fprintf(stdout,"Populating training set on proc %i...\n",rank);

    gsl_vector_complex *wv;
    double *params;
    int start_ind, end_ind, global_i, matrix_size;


    params = new double[ts.param_dim]; // for TF2 gravitational wave model this is (mass 1, mass 2)
    wv     = gsl_vector_complex_alloc(xQuad->size);

    // -- decide which chunk of TS to compute on this proc -- //
    if(ts.distributed){
        start_ind   = ts.mystart[rank];
        end_ind     = ts.myend[rank];
        matrix_size = ts.matrix_sub_size[rank];
        fprintf(stdout,"start ind is %i and end ind is %i\n",start_ind,end_ind);
    }
    else{
        start_ind   = 0;
        matrix_size = ts.ts_size;
        end_ind     = ts.ts_size;
    }


    // *** BEGIN MODEL SPECIFIC SECTION *** //
    // This is where a new model should go...add to the list and loop over paramters //
    if(strcmp(ts.model,"TaylorF2_PN3pt5") == 0){
        fprintf(stdout,"Using the TaylorF2 spa approximant to PN=3.5\n");

        for(int i = 0; i < matrix_size; i++){

            global_i = start_ind + i;

            // TODO: this should be handled in ts class whenever thats implimented
            for(int j = 0; j < ts.param_dim; j++){
                params[j] = ts.params[global_i][j] * ts.param_scale[j];
            }

            TF2_FullWaveform(wv,params,xQuad,1.0,3.5); // amp = 1.0 and PN order 3.5
            gsl_matrix_complex_set_row(TS_gsl,i,wv);
        }

    }
    else{
        std::cerr << "Approximant not supported!" << std::endl;
        exit(1);
    }    
    // *** END MODEL SPECIFIC SECTION *** //



    // -- Normalize training space here -- //
    fprintf(stdout,"Normalizing training set...\n");
    NormalizeMatrixRows(TS_gsl,wQuad);

    delete[] params;
    gsl_vector_complex_free(wv);
}

void MGS(gsl_vector_complex *ru, gsl_vector_complex *ortho_basis,const gsl_matrix_complex *RB_space,const gsl_vector_complex *wQuad, const int dim_RB)
{
/*  Modified GS routine. 

   Input:  RB_space:    an existing orthonormal set of basis vectors
           ortho_basis: basis we shall orthogonalize agains RB_space
           wQuad:       quadrature weights for inner products
           dim_RB:      number of elements currently in RB_space
   Output: ortho_basis: orthonormal basis vector
           ru:           1-by-(dim_RB+1) slice of the R matrix (from QR = A)
*/
          

    int quad_num = RB_space->size2;
    gsl_complex L2_proj, tmp;
    gsl_vector_complex *basis;


    basis = gsl_vector_complex_alloc(quad_num);

    gsl_vector_complex_set_zero(ru); // if not done, R matrix fills up below diagonal with .5 instead of 0


    for(int i = 0; i < dim_RB; i++)
    {
        gsl_matrix_complex_get_row(basis,RB_space,i);

        /* --- ortho_basis = ortho_basis - L2_proj*basis; --- */
        L2_proj = WeightedInner(wQuad,basis,ortho_basis);
        gsl_vector_complex_set(ru,i,L2_proj);
        gsl_vector_complex_scale(basis,L2_proj); // basis <- basis*L2_proj
        gsl_vector_complex_sub(ortho_basis,basis); // ortho_basis <- ortho_basis - basis

    }

    double nrm = GetNorm_double(ortho_basis,wQuad);
    gsl_complex nrmc;
    GSL_SET_COMPLEX(&nrmc,nrm,0.0);
    gsl_vector_complex_set(ru,dim_RB,nrmc);

    NormalizeVector(ortho_basis,wQuad);

    gsl_vector_complex_free(basis);

}


void IMGS(gsl_vector_complex *ru, gsl_vector_complex *ortho_basis,const gsl_matrix_complex *RB_space,const gsl_vector_complex *wQuad, const int dim_RB)
{
/*  Iterated modified GS routine. 

   Input:  RB_space:    an existing orthonormal set of basis vectors
           ortho_basis: basis we shall orthogonalize agains RB_space
           wQuad:       quadrature weights for inner products
           dim_RB:      number of elements currently in RB_space (ortho_basis is dim_RB+1 element)
   Output: ortho_basis: orthonormal basis vector
           ru:           1-by-(dim_RB+1) slice of the R matrix (from QR = A)

   Hoffmann, "ITERATIVE ALGORITHMS FOR GRAM-SCHMIDT ORTHOGONALIZATION"

*/

    double ortho_condition = .5; // IMGS stopping condition (HARD CODED!!)

    int quad_num     = RB_space->size2;
    int r_size       = ru->size;
    double nrm_prev  = GetNorm_double(ortho_basis,wQuad);
    bool flag        = false;
    int iter         = 0;
    gsl_vector_complex *e,*r_last;
    double nrm_current;
    gsl_complex nrmc_current;

    // --- allocate memory --- //
    e = gsl_vector_complex_alloc(quad_num);
    r_last = gsl_vector_complex_alloc(r_size);

    gsl_vector_complex_memcpy(e,ortho_basis);
    NormalizeVector(e,wQuad);
    gsl_vector_complex_set_zero(ru); // if not done, R matrix fills up below diagonal with .5 instead of 0

    // TODO: code below looks to have some redundent parts -- after working cleanup

    while(!flag)
    {
        gsl_vector_complex_memcpy(ortho_basis,e);
        gsl_vector_complex_set_zero(r_last);

        MGS(r_last,ortho_basis,RB_space,wQuad,dim_RB);

        gsl_vector_complex_add(ru,r_last);
        nrmc_current = gsl_vector_complex_get(r_last,dim_RB);
        nrm_current = GSL_REAL(nrmc_current);

        gsl_vector_complex_scale(ortho_basis,nrmc_current);


        if( nrm_current/nrm_prev <= ortho_condition )
        {
            nrm_prev = nrm_current;
            iter = iter + 1;
            gsl_vector_complex_memcpy(e,ortho_basis);

        }
        else
        {
            flag = true;
        }

        nrm_current  = GetNorm_double(ortho_basis,wQuad);
        GSL_SET_COMPLEX(&nrmc_current,nrm_current,0.0);
        gsl_vector_complex_set(ru,dim_RB,nrmc_current);

        NormalizeVector(ortho_basis,wQuad);

    }

    gsl_vector_complex_free(e);
    gsl_vector_complex_free(r_last);

}

int FindRowIndxRank(const int global_row_indx,const TrainSet ts)
{
    int row_rank = 0;
    bool found_rank = false;

    if(ts.distributed){

        while(found_rank == false)
        {
            if( global_row_indx >= ts.mystart[row_rank] && global_row_indx <= ts.myend[row_rank] ){
                found_rank = true;
            }
            else{
                row_rank = row_rank + 1;
            }
        }
    }

    return row_rank;

}

void WriteGreedyInfo(const int dim_RB, const gsl_matrix_complex *RB_space, const gsl_matrix_complex *R_matrix, const double *app_err, const int *sel_rows, const TrainSet ts, const char * output_dir,const char *datatype)
{
    FILE *rb_real_data, *rb_imag_data, *r_real_data, *r_imag_data, *err_data, *pts_data;
    FILE *rb_data, *r_data;
    char rb_real_filename[100];
    char rb_imag_filename[100];
    char r_real_filename[100];
    char r_imag_filename[100];
    char err_filename[100];
    char pts_filename[100];
    char rb_filename[100];
    char r_filename[100];

    if(strcmp(datatype,"txt") == 0){
        strcpy(rb_real_filename,output_dir);
        strcat(rb_real_filename,"/Basis_real.txt");
        strcpy(rb_imag_filename,output_dir);
        strcat(rb_imag_filename,"/Basis_imag.txt");
        strcpy(r_real_filename,output_dir);
        strcat(r_real_filename,"/R_real.txt");
        strcpy(r_imag_filename,output_dir);
        strcat(r_imag_filename,"/R_imag.txt");
    } 
    else if(strcmp(datatype,"bin") == 0){
        strcpy(rb_filename,output_dir);
        strcat(rb_filename,"/Basis.bin");
        strcpy(r_filename,output_dir);
        strcat(r_filename,"/R.bin");
    }
    else{
        fprintf(stderr,"file type not supported");
        exit(1);
    }

    strcpy(err_filename,output_dir);
    strcat(err_filename,"/ApproxErrors.txt");
    strcpy(pts_filename,output_dir);
    strcat(pts_filename,"/GreedyPoints.txt");

    //--- write errors and greedy points to text file ---//
    err_data = fopen(err_filename,"w");
    pts_data = fopen(pts_filename,"w");
    for(int i = 0; i < dim_RB ; i++)
    {
        fprintf(err_data,"%1.14f\n",app_err[i]);
        fprintf(pts_data,"%1.14f %1.14f\n",ts.params[sel_rows[i]][0],ts.params[sel_rows[i]][1]);
    }
    fclose(err_data);
    fclose(pts_data);

    //--- write R and RB to file ---//
    if(strcmp(datatype,"txt") == 0){
        rb_real_data = fopen(rb_real_filename,"w");
        gsl_matrix_complex_fprintf_part(rb_real_data,RB_space,"real");
        fclose(rb_real_data);
        rb_imag_data = fopen(rb_imag_filename,"w");
        gsl_matrix_complex_fprintf_part(rb_imag_data,RB_space,"imag");
        fclose(rb_imag_data);
        r_real_data = fopen(r_real_filename,"w");
        gsl_matrix_complex_fprintf_part(r_real_data,R_matrix,"real");
        fclose(r_real_data);
        r_imag_data = fopen(r_imag_filename,"w");
        gsl_matrix_complex_fprintf_part(r_imag_data,R_matrix,"imag");
        fclose(r_imag_data);
    }
    else{
        rb_data = fopen(rb_filename,"w");
        gsl_matrix_complex_fwrite(rb_data,RB_space);
        fclose(rb_data);
        r_data = fopen(r_filename,"w");
        gsl_matrix_complex_fwrite(r_data,R_matrix);
        fclose(r_data);
    }

}

void WriteWaveform(const double *xQuad,const gsl_matrix_complex *TS_gsl,const int indx)
{

    FILE *data;
    char filename[] = "MockWaveform.txt";
    data = fopen(filename,"w");
    for(int cols = 0; cols < TS_gsl->size2 ; cols++)
    {
        fprintf(data,"%1.12e %1.12e %1.12e\n",xQuad[cols], GSL_REAL(gsl_matrix_complex_get(TS_gsl,indx,cols)),GSL_IMAG( gsl_matrix_complex_get(TS_gsl,indx,cols) ));
    }
    fclose(data);

}


void GreedyWorker(const int rank, const int max_RB,const int seed_global, const gsl_vector_complex *wQuad,const double tol, const gsl_matrix_complex *A, const TrainSet ts)
{

// worker routine for computing computationally intensive part of greedy //

    int continue_work = 1;
    int dim_RB = 1;
    double cols = wQuad->size;
    int worst_global, worst_local, worst_rank;
    double rb_inc_r, rb_inc_i, tmp;
    double *errors, *vec_real, *vec_imag;
    gsl_complex tmpc;
    gsl_vector_complex *last_rb, *row_vec;
    gsl_matrix_complex *project_coeff;

    last_rb       = gsl_vector_complex_alloc(cols);
    row_vec        = gsl_vector_complex_alloc(cols);
    project_coeff = gsl_matrix_complex_alloc(max_RB,ts.matrix_sub_size[rank]);
    errors        = (double *)malloc(ts.matrix_sub_size[rank]*sizeof(double));
    vec_real      = (double *)malloc(cols*sizeof(double));
    vec_imag      = (double *)malloc(cols*sizeof(double));

    int *worst_workers_mpi = NULL;
    double *worst_errs_mpi = NULL;

    fprintf(stdout,"I'm worker %i and I was given %i matrix elements from %i to %i\n",rank,ts.matrix_sub_size[rank],ts.mystart[rank],ts.myend[rank]-1);

    // -- pass seed back to master -- //
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&worst_rank, 1, MPI_INT,0,MPI_COMM_WORLD);
    // -- return seed to master -- //
    if( (worst_rank-1) == rank){
        worst_local = seed_global - ts.mystart[rank];
        gsl_matrix_complex_get_row(row_vec,A,worst_local);
        gsl_vector_complex_parts(vec_real,vec_imag,row_vec);
        MPI_Send(vec_real,cols, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        MPI_Send(vec_imag,cols, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }

    while(continue_work == 1)
    {
        // -- wait for new basis and start next sweep -- //
        MPI_Barrier(MPI_COMM_WORLD);

        // -- receive new rb -- //
        MPI_Bcast(vec_real, cols, MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Bcast(vec_imag, cols, MPI_DOUBLE,0,MPI_COMM_WORLD);
        make_gsl_vector_complex_parts(vec_real,vec_imag,last_rb);

        // Compute overlaps of pieces of A with rb_new //
        for(int i = 0; i < ts.matrix_sub_size[rank]; i++)
        {
            gsl_matrix_complex_get_row(row_vec,A,i);
            tmpc = WeightedInner(wQuad,last_rb,row_vec);
            gsl_matrix_complex_set(project_coeff,dim_RB-1,i,tmpc);

            tmp = 0;
            for(int j = 0; j < dim_RB; j++)
            {
                tmpc = gsl_matrix_complex_get(project_coeff,j,i);
                tmp = tmp + gsl_complex_abs(tmpc)*gsl_complex_abs(tmpc);
            }
            errors[i] = 1.0 - tmp;
        }

        // -- find worst error here -- //
        tmp = 0.0;
        for(int i = 0; i < ts.matrix_sub_size[rank]; i++)
        {
            if(tmp < errors[i])
            {
                tmp = errors[i];
                worst_local = i;                // this will retun matrix local (worker's) row index
                worst_global = ts.mystart[rank] + i; // this will return matrix's global row index 
            }
        }

        // -- pass worst error and index to master --//
        MPI_Gather(&worst_global,1,MPI_INT,worst_workers_mpi,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Gather(&tmp,1,MPI_DOUBLE,worst_errs_mpi,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

        // -- recieve info about which worker proc has next basis -- //
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(&worst_rank, 1, MPI_INT,0,MPI_COMM_WORLD);
 
        // -- return basis to master -- //
        if( (worst_rank-1) == rank){
            gsl_matrix_complex_get_row(row_vec,A,worst_local);
            gsl_vector_complex_parts(vec_real,vec_imag,row_vec);
            MPI_Send(vec_real,cols, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Send(vec_imag,cols, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(&continue_work, 1, MPI_INT,0,MPI_COMM_WORLD);
        dim_RB = dim_RB + 1;
    }

    free(worst_workers_mpi); // free NULL pointers?
    free(worst_errs_mpi);

    gsl_vector_complex_free(last_rb);
    gsl_vector_complex_free(row_vec);
    gsl_matrix_complex_free(project_coeff);
    free(vec_real);
    free(vec_imag);
    free(errors);
}

void GreedyMaster(const int size, const int max_RB, const int seed,const gsl_vector_complex *wQuad,const double tol, const TrainSet ts, const char * output_dir, const char *output_data_format)
{
// Input: 
//          A: gsl matrix of solutions (each row is a solutions, cols are quadrature samples)
//          seed: first greedy pick (global indexing)
//          tol: approximation tolerance
//          size: number of procs. if 1 serial mode assumed
//
//   Output:
//          sel_rows: row index defining reduced basis. sel_rows[0] = seed
//          dim_RB: number of greedy_points
//

    fprintf(stdout,"Starting greedy algorithm...\n");

    const int rows = ts.ts_size;  // number of rows to approximate
    const int cols = wQuad->size; // samples (for quadrature)
    int *greedy_points;           // selectted greedy points (row selection)
    double *greedy_err;           // approximate error
    clock_t start, end;           // for algorithm timing experiments
    double alg_time;              // for algorithm timing experiments
    double worst_err;             // errors in greedy sweep
    int worst_app, worst_worker, worst_rank;             // worst error stored
    double rb_inc_r, rb_inc_i;
    double *vec_real, *vec_imag;
    int *worst_workers_mpi;
    double *worst_errs_mpi;
    gsl_complex tmpc;
    int dummy_mpi_int        = -1;
    double dummy_mpi_double  = -1.0;
    int continue_work = 1;
    int dim_RB       = 1;

    gsl_vector_complex *ortho_basis, *ru;
    gsl_matrix_complex *RB_space, *R_matrix;

    // --- this memory should be freed here --- //
    ortho_basis       = gsl_vector_complex_alloc(cols);
    vec_real          = (double *)malloc(cols*sizeof(double));
    vec_imag          = (double *)malloc(cols*sizeof(double));
    worst_workers_mpi = (int *)malloc(size*sizeof(int));
    worst_errs_mpi    = (double *)malloc(size*sizeof(double));
    RB_space          = gsl_matrix_complex_alloc(max_RB,cols); 
    R_matrix          = gsl_matrix_complex_alloc(max_RB,max_RB);
    greedy_points     = (int *)malloc(max_RB*sizeof(int));
    greedy_err        = (double *)malloc(max_RB*sizeof(double));
    ru                = gsl_vector_complex_alloc(max_RB);

    // --- initialize algorithm with seed --- //
    int seed_rank = FindRowIndxRank(seed,ts);
    fprintf(stdout,"seed index %i on proc rank %i\n",seed,seed_rank);

    // -- request seed from worker -- //
    MPI_Barrier(MPI_COMM_WORLD);
    seed_rank = seed_rank+1;
    MPI_Bcast(&seed_rank, 1, MPI_INT,0,MPI_COMM_WORLD);
    MPI_Recv(vec_real, cols, MPI_DOUBLE, seed_rank, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    MPI_Recv(vec_imag, cols, MPI_DOUBLE, seed_rank, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    make_gsl_vector_complex_parts(vec_real,vec_imag,ortho_basis);
    gsl_matrix_complex_set_row(RB_space,0,ortho_basis);


    GSL_SET_COMPLEX(&tmpc,1.0,0.0); // assumes normalized solutions
    gsl_matrix_complex_set(R_matrix,0,0,tmpc);
    greedy_points[0] = seed;
    dim_RB           = 1;
    greedy_err[0]    = 1.0;

    // --- Continue approximation until tolerance satisfied --- //
    start = clock();
    while(continue_work == 1)
    {
        gsl_vector_complex_parts(vec_real,vec_imag,ortho_basis);

        // -- send last orthonormal rb to work procs -- //
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(vec_real, cols, MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Bcast(vec_imag, cols, MPI_DOUBLE,0,MPI_COMM_WORLD);

        // -- gather worst (local) info from workers -- //
        // TODO: worst errors can be passed along with basis, no need to gather
        MPI_Gather(&dummy_mpi_int,1,MPI_INT,worst_workers_mpi,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Gather(&dummy_mpi_double,1,MPI_DOUBLE,worst_errs_mpi,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

        // -- find worst rb amongst all workers -- //
        worst_err = 0.0;
        for(int i = 0; i < size - 1; i++){
            if(worst_err < worst_errs_mpi[i+1])
            {
                worst_err  = worst_errs_mpi[i+1];
                worst_app  = worst_workers_mpi[i+1];
                worst_rank = i+1;
            }
        }

        // -- tell all workers which one has largest error -- //
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(&worst_rank, 1, MPI_INT,0,MPI_COMM_WORLD);

        // -- receive row basis from worker proc worst_rank -- //
        MPI_Recv(vec_real, cols, MPI_DOUBLE, worst_rank, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MPI_Recv(vec_imag, cols, MPI_DOUBLE, worst_rank, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        make_gsl_vector_complex_parts(vec_real,vec_imag,ortho_basis);

        // -- decide if another greedy sweep is needed, alert workers -- //
        if( (dim_RB+1 == max_RB) || worst_err < tol){
            continue_work = 0;
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(&continue_work, 1, MPI_INT,0,MPI_COMM_WORLD);

        // --- record worst approximated row element (basis) --- //
        greedy_points[dim_RB] = worst_app;
        greedy_err[dim_RB] = worst_err;

        // --- add worst approximated solution/row to basis set --- //
        IMGS(ru,ortho_basis,RB_space,wQuad,dim_RB); // IMGS SHOULD BE DEFAULT
        //MGS(ru,ortho_basis,RB_space,wQuad,dim_RB);
        gsl_matrix_complex_set_row(R_matrix,dim_RB,ru);
        gsl_matrix_complex_set_row(RB_space,dim_RB,ortho_basis);
        dim_RB = dim_RB + 1;

        fprintf(stdout,"RB dimension %i || Current row selection %i || Approximation error %1.14e\n",dim_RB,worst_app,worst_err);

    }
    end = clock();

    alg_time = ((double) (end - start)/CLOCKS_PER_SEC);
    fprintf(stdout,"Building approximation space took %f seconds\n",alg_time);
    dim_RB = dim_RB - 1;

    // -- output relevant information -- //
    WriteGreedyInfo(dim_RB,RB_space,R_matrix,greedy_err,greedy_points,ts,output_dir,output_data_format);


    gsl_vector_complex_free(ortho_basis);
    free(vec_real);
    free(vec_imag);
    free(worst_workers_mpi);
    free(worst_errs_mpi);
    gsl_matrix_complex_free(RB_space); 
    gsl_matrix_complex_free(R_matrix);
    free(greedy_points);
    free(greedy_err);
    gsl_vector_complex_free(ru);

}

void Greedy(const int seed,const int max_RB, const gsl_matrix_complex *A,const gsl_vector_complex *wQuad,const double tol,const TrainSet ts, const char * output_dir, const char *output_data_format)
{
// Input: 
//          A: gsl matrix of solutions (each row is a solution, cols are quadrature samples)
//          seed: first greedy pick
//          tol: approximation tolerance
//
// Output:
//          sel_rows: row index defining reduced basis. sel_rows[0] = seed
//          dim_RB: number of greedy_points
//

    fprintf(stdout,"Starting greedy algorithm in serial mode...\n");

    const int rows = A->size1; // number of rows to approximate
    const int cols = A->size2; // samples (for quadrature)
    int *greedy_points;        // selectted greedy points (row selection)
    double *greedy_err;        // approximate error
    clock_t start, end;        // for algorithm timing experiments
    double alg_time;           // for algorithm timing experiments
    double tmp,worst_err;      // errors in greedy sweep
    int worst_app;             // worst error stored
    gsl_complex tmpc;          // worst error temp
    bool continue_work = true;


    gsl_vector_complex *ts_el, *last_rb, *ortho_basis, *ru;
    gsl_matrix_complex *RB_space, *R_matrix;
    double *errors;                       // approximation errors with RB space of dimension dim_RB
    gsl_matrix_complex *project_coeff; // h = coeff_i e_i is approximation we seek


    // --- this memory should be freed here --- //
    ts_el         = gsl_vector_complex_alloc(cols);
    last_rb       = gsl_vector_complex_alloc(cols);
    ortho_basis   = gsl_vector_complex_alloc(cols);
    ru = gsl_vector_complex_alloc(max_RB);
    errors        = (double *)malloc(rows*sizeof(double));
    project_coeff = gsl_matrix_complex_alloc(max_RB,rows);
    greedy_points = (int *)malloc(max_RB*sizeof(int));
    greedy_err    = (double *)malloc(max_RB*sizeof(double));
    RB_space      = gsl_matrix_complex_alloc(max_RB,cols); 
    R_matrix      = gsl_matrix_complex_alloc(max_RB,max_RB);

    // --- initialize algorithm with seed --- //
    gsl_matrix_complex_get_row(ts_el,A,seed);
    gsl_matrix_complex_set_row(RB_space,0,ts_el);
    GSL_SET_COMPLEX(&tmpc,1.0,0.0);
    gsl_matrix_complex_set(R_matrix,0,0,tmpc);  // assumes normalized solutions

    greedy_points[0] = seed;
    int dim_RB       = 1;
    greedy_err[0]    = 1.0;

    // --- Continue approximation until tolerance satisfied --- //
    start = clock();
    while(continue_work)
    {

        gsl_matrix_complex_get_row(last_rb,RB_space,dim_RB-1); // get last computed basis
        worst_err = 0.0;

        // --- Loop over training set ---//
        for(int i = 0; i < rows; i++)
        {

            gsl_matrix_complex_get_row(ts_el,A,i);
            tmpc = WeightedInner(wQuad,last_rb,ts_el);
            gsl_matrix_complex_set(project_coeff,dim_RB-1,i,tmpc);

            tmp = 0;
            for(int j = 0; j < dim_RB; j++)
            {
               tmpc = gsl_matrix_complex_get(project_coeff,j,i);
               tmp = tmp + gsl_complex_abs(tmpc)*gsl_complex_abs(tmpc);
            }

            errors[i] = 1.0 - tmp;

            if(worst_err < errors[i])
            {
                worst_err = errors[i];
                worst_app = i;
            }

        }

        // --- add worst approximated element to basis --- //
        greedy_points[dim_RB] = worst_app;
        greedy_err[dim_RB] = worst_err;

        // -- decide if another greedy sweep is needed -- //
        if( (dim_RB+1 == max_RB) || (worst_err < tol) || (ts.ts_size == dim_RB) ){
            continue_work = false;
        }


        // --- add worst approximated solution to basis set --- //
        gsl_matrix_complex_get_row(ortho_basis,A,worst_app);

        IMGS(ru,ortho_basis,RB_space,wQuad,dim_RB); // IMGS SHOULD BE DEFAULT

        //MGS(ru,ortho_basis,RB_space,wQuad,dim_RB);
        gsl_matrix_complex_set_row(RB_space,dim_RB,ortho_basis);
        gsl_matrix_complex_set_row(R_matrix,dim_RB,ru);
        dim_RB = dim_RB + 1;


        fprintf(stdout,"RB dimension %i || Current row selection %i || Approximation error %1.14e\n",dim_RB,worst_app,worst_err);

    }
    end = clock();

    alg_time = ((double) (end - start)/CLOCKS_PER_SEC);
    fprintf(stdout,"Building approximation space took %f seconds\n",alg_time);
    dim_RB = dim_RB - 1;

    // -- output relevant information -- //
    WriteGreedyInfo(dim_RB,RB_space,R_matrix,greedy_err,greedy_points,ts,output_dir,output_data_format);


    gsl_vector_complex_free(ts_el);
    gsl_vector_complex_free(last_rb);
    gsl_vector_complex_free(ortho_basis);
    gsl_vector_complex_free(ru);
    free(errors);
    gsl_matrix_complex_free(project_coeff);
    free(greedy_points);
    free(greedy_err);
    gsl_matrix_complex_free(RB_space); // this and R_matrix should be written to file
    gsl_matrix_complex_free(R_matrix);

}

static int fcount_pts(const char *infile) {
    //const char *c_str = infile.c_str();
    //FILE *fvec = fopen(c_str, "r");
    std::ifstream fvec(infile);
    // counts the number of lines in the given file not beginning with "#"
    // to get the number of points
    std::string line;
    int points = 0;
    while (std::getline(fvec, line))
    {
        if (line.compare(0, 1, "#") != 0)
            ++points;

    }
    //fclose(fvec);

    return points;
}

int main (int argc, char **argv) {

    // --- setup MPI info ---//
    MPI::Init(argc, argv);

    int rank = 0;  // needed for serial mode too
    int size = 1;  // needed for serial mode too

    // Get number of procs this job is using (size) and unique rank of the processor this thread is running on //
    // Ex: if executed with "mpirun -np 2", size = 2. "CPU1" will be 0 and "CPU2" will be 1 //
    rank = MPI::COMM_WORLD.Get_rank();
    size = MPI::COMM_WORLD.Get_size();

    char name[MPI_MAX_PROCESSOR_NAME];
    int len;
    memset(name,0,MPI_MAX_PROCESSOR_NAME);
    MPI::Get_processor_name(name,len);
    memset(name+len,0,MPI_MAX_PROCESSOR_NAME-len);

    std::cout << "Number of tasks = " << size << " My rank = " << rank << " My name = " << name << "." << std::endl;

    //----- Checking the number of Variables passed to the Executable -----//
    if (argc != 2) {
        std::cerr << "Arguments: 1. location of a cfg configuration/parameter file (ends in .cfg)" << std::endl;
        exit(0);
    }
    std::cout << "parameter file is: " << argv[1] << std::endl;


    //--- Read input (config) file. If there is an error, report and exit ---//
    // TODO: with MPI, multiple procs reading same paramter file... seems bad //
    libconfig::Config cfg;
    try{
      cfg.readFile(argv[1]);
    }
    catch(const libconfig::FileIOException &fioex){
      std::cerr << "I/O error while reading file." << std::endl;
      return(EXIT_FAILURE);
    }
    catch(const libconfig::ParseException &pex){
      std::cerr << "Parse error " << std::endl;
      return(EXIT_FAILURE);
    }

    // Variables which need to be decarled, but not set in parameter file //
    gsl_matrix_complex *TS_gsl;
    gsl_vector_complex *wQuad;
    gsl_vector *xQuad;
    TrainSet ts;
    char shell_command[100];
    int gsl_status;
    bool cfg_status;
    int ts_size;

    // run settings - these MUST be set in the parameter file //
    bool load_from_file = cfg.lookup("load_from_file");  // load training points from file instead (file name used is below)
    const int quad_type = cfg.lookup("quad_type");       // 0 = Gaussian quadrature, 1 = Riemann, 2 = user-defined (quadrature file required)
    const int seed      = cfg.lookup("seed");            // greedy algorithm seed
    const double tol    = cfg.lookup("tol");             // greedy algorithm tolerance ( \| app \|^2)
    int max_RB          = cfg.lookup("max_RB");          // estimated number of RB (reasonable upper bound)
    bool weighted_inner = cfg.lookup("weighted_inner");  // whether or not the inner product to use includes a weight W(x): \int W(x) f(x) g(x)
    int param_dim       = cfg.lookup("param_dim");       // number of paramteric dimensions (currently supports 2)
    const char * model_name;                             // mame of model -- used to select the appropriate model in FillTrainingSet
    const char * output_dir;                             // folder to put all output files
    const char * output_data_format;                     // format of output files (text or gsl binary supported)

    // TODO: should be function in TS routine for ND parameter spaces
    ts.param_scale = (double *)malloc(param_dim*sizeof(double)); 
    char scale_str[20];
    for(int i = 0; i < param_dim; i++){
        snprintf(scale_str, 20, "p%d_scale", i+1);
        ts.param_scale[i] =  cfg.lookup(scale_str);   // scale each params[j][i] so that model evaluated at param_scale[i] * params[j][i]      
    }


    // run settings - MAY need to be set for SetupQuadratureRule //
    // x_min;              // lower value x_min (physical domain) --> needed if quad_type != 2
    // x_max;              // upper value x_max (physical domain) --> needed if quad_type != 2
    // quad_points;        // total number of quadrature points   --> needed if quad_type != 2
    // fvec_file_name;     // file name for vector of quadrature points  --> needed if quad_type = 2
    // weight_file_name;   // file name of weights --> needed if weighted_inner = true

    // run settings - MAY need to be set //
    const char *ts_file_name;          // this is required if load_from_file = true
    double *params_low, *params_high;  // this is required if load_from_file = false. lower/upper interval of each parameter 
    int *params_num;                   // this is required if load_from_file = false. Number of samplings in the interval [param_low,param_high]

    // lookup additional required run parameters from configuration file //
    if(cfg.lookupValue("model_name",model_name) 
       && cfg.lookupValue("output_dir",output_dir) 
       && cfg.lookupValue("output_data_format",output_data_format)){
        fprintf(stdout,"Successfully loaded model name, output directory location and output data format type\n");
    }
    else{
        fprintf(stderr,"Failed to load either model name, output directory location or output data format type\n");
        exit(1);
    }


    // returns wQuad and xQuad, configuration file = argv[1]  //
    // note: this returns the full, not reduced, quadrature rule //
    SetupQuadratureRule(&wQuad,&xQuad,quad_type,weighted_inner,argv[1]);

    const libconfig::Setting& root = cfg.getRoot();

    // -- build training set (builds list of paramter values in ts.params) -- //
    if (load_from_file){
        cfg_status = cfg.lookupValue("ts_file", ts_file_name);
        if (!cfg_status){
            fprintf(stderr, "ts_file not found in config file\n");
            exit(1);
        }
        ts_size = fcount_pts(ts_file_name); // need to get training set size from file (number of rows)
        std::cout << "training set file found to " << ts_size << " parameter samples" << std::endl;
    }
    else{
        params_num       = (int *)malloc(param_dim*sizeof(int));
        params_low       = (double *)malloc(param_dim*sizeof(double));
        params_high      = (double *)malloc(param_dim*sizeof(double));

        // note about libconfig: if destination data type does not match expected on from configuration file an error is thrown //

        // read in params_num and determine ts_size //
        ts_size = 1;
        libconfig::Setting& params_num_s = root["params_num"];
        int count_array = params_num_s.getLength();
        if(count_array != param_dim){
            fprintf(stderr,"elements in params_num defined in configuration file does equal param_dim\n");
            exit(1);
        }
        for(int i = 0; i < param_dim; i++){
            params_num[i] = params_num_s[i];
            ts_size       = params_num[i]*ts_size; // assumes tensor product structure on training set
        }

        // read in params_low //
        libconfig::Setting& params_low_s = root["params_low"];
        count_array = params_low_s.getLength();
        if(count_array != param_dim){
            fprintf(stderr,"elements in params_low defined in configuration file does equal param_dim\n");
            exit(1);
        }
        for(int i = 0; i < param_dim; i++){
            params_low[i] = params_low_s[i];
        }

        // read in params_low //
        libconfig::Setting& params_high_s = root["params_high"];
        count_array = params_high_s.getLength();
        if(count_array != param_dim){
            fprintf(stderr,"elements in params_high defined in configuration file does equal param_dim\n");
            exit(1);
        }
        for(int i = 0; i < param_dim; i++){
            params_high[i] = params_high_s[i];
        }
    }

    // -- Finished reading from configuration file. Start algorithm... //



    // Creating Run Directory //
    if(size == 1 || rank == 0){

        strcpy(shell_command, "mkdir -p -m700 ");
        strcat(shell_command, output_dir);
        system(shell_command);

        snprintf(shell_command,100,"cp %s %s",argv[1],output_dir);
        system(shell_command);
    }

    // this function sets members of ts structure and allocates memory for set of parameters //
    ts_alloc(ts_size, param_dim, model_name, ts);

    // returns ts.params filled with values //
    if(load_from_file){
        BuildTS_FromFile(ts_file_name,ts);
        //BuildTS_FromFile_ND(ts_file_name,ts); // N-dimensional routine in testing... should work
    }
    else{
        BuildTS_TensorProduct(params_num,params_low,params_high,ts);
    }

    // Build training space by evaluating model at ts.params. Then run the greedy algorithm //
    if(size == 1) // only 1 proc requested (serial mode)
    {
        TS_gsl = gsl_matrix_complex_alloc(ts.ts_size,xQuad->size); // GSL error handler will abort if too much requested
        FillTrainingSet(TS_gsl,xQuad,wQuad,ts,rank);               // size=1  => rank=0. 5th argument is the rank
        Greedy(seed,max_RB,TS_gsl,wQuad,tol,ts,output_dir,output_data_format);
    }
    else
    {
        // -- split matrix TS_gsl among worker nodes. Assumes for-loop is "<" for this choice of myend -- //
        SplitTrainingSet(size,ts);

        if(rank != 0){
            TS_gsl = gsl_matrix_complex_alloc(ts.matrix_sub_size[rank-1],xQuad->size);
            FillTrainingSet(TS_gsl,xQuad,wQuad,ts,rank-1);
        }

        fprintf(stdout,"Finished distribution of training set\n");

        if(rank == 0){
            GreedyMaster(size,max_RB,seed,wQuad,tol,ts,output_dir,output_data_format);
        }
        else{
            GreedyWorker(rank-1,max_RB,seed,wQuad,tol,TS_gsl,ts);
            gsl_matrix_complex_free(TS_gsl);
        }

    }


    if(rank == 0)
    {

        // -- output quadrature weights -- //
        FILE *outfile;

        strcpy(shell_command, output_dir);
        strcat(shell_command,"/quad_weights.txt");

        outfile = fopen(shell_command,"w");

        for(int i = 0; i < wQuad->size ; i++) {
            fprintf(outfile,"%1.14f\n",GSL_REAL(gsl_vector_complex_get(wQuad,i)));
        }
        fclose(outfile);

        // -- output some waveform for diagnostics -- //
        //if(size == 1){
        //    WriteWaveform(xQuad->data,TS_gsl,0); // for comparison with other codes
        //    gsl_matrix_complex_free(TS_gsl);
        //}
        //WriteTrainingSet(ts);
    }


    gsl_vector_complex_free(wQuad);
    gsl_vector_free(xQuad);
    free(ts.param_scale);

    // free params matrix. TODO: check with valgrind that this is correct by looking for memory leaks //
    for(int j = 0; j < ts.ts_size; j++){ 
      free(ts.params[j]);
    }
    free(ts.params);

    if(ts.distributed){
        free(ts.mystart);
        free(ts.myend);
        free(ts.matrix_sub_size);
    }

    if (load_from_file == false){
        free(params_num);
        free(params_low);
        free(params_high);
    }

    // Tell the MPI library to release all resources it is using
    MPI::Finalize();

}

