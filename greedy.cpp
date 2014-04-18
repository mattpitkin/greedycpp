// AUTHOR :  Scott Field 
//           sfield@umd.edu
//
// DATE: Jan 20, 2013
//
// PURPOSE: mpi version of greedy (pivoted MGS) algorithm 


//#include "nr3.h"
#include "gauss_wgts.h"

#include <mpi.h> // MPI STUFF

//#include <boost/numeric/ublas/vector.hpp>
//#include <boost/numeric/ublas/io.hpp>
#include <iostream>
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
#include "psd.h"
#include "TrainingSet.h"

/*-- constants --*/
double c = 299792458.0;                             // speed of light
double c2 = c*c;
double mass_to_sec = 1.98892e30*6.67428E-11/(c2*c); // Conversion factor to change Solar Mass units to seconds. assumes G = 6.67428E-11 and SolarMass = 1.98892e30
double DL = 3.08568025E22;                          // Luminosity distance in meters, set to one megaparsec


/* --- fill array with linear spacing --- */
void Linspace(const int &n, const double &a, const double &b, double *SomeArray)
{
    double factor = (b-a)/(double)(n-1);
    for(int i=0;i<n;i++)
    {
        SomeArray[i] = a + (double)i*factor;
    }
}

void gsl_vector_complex_parts(double *v_real, double *v_imag, const gsl_vector_complex * v_gsl)
{

    for(int i = 0; i < v_gsl->size; i++)
    {
        v_real[i] = GSL_REAL(gsl_vector_complex_get(v_gsl,i));
        v_imag[i] = GSL_IMAG(gsl_vector_complex_get(v_gsl,i));
    }
}

void make_gsl_vector_complex_parts(const double *v_real, const double *v_imag,gsl_vector_complex * v_gsl)
{
    gsl_complex tmpc;

    for(int i = 0; i < v_gsl->size; i++){
        GSL_SET_COMPLEX(&tmpc,v_real[i],v_imag[i]);
        gsl_vector_complex_set(v_gsl,i,tmpc);
    }
}

void ReimannQuad(const double a,const double b,double *xQuad,double * wQuad,const int freq_points)
{

    Linspace(freq_points, a, b, xQuad);

    for(int i = 0; i < freq_points; i++)
    {
        wQuad[i] = (b-a)/( (double) (freq_points-1) );
    }

}

void OutputArray(const int n, double *list)
{
    for(int i=0;i<n;i++){
        std::cout << list[i] << std::endl;
    }
}

void MakeQuadratureRule(gsl_vector_complex *wQuad_c, gsl_vector *xQuad_c, const double a, const double b, const int freq_points,const int quad_type)
{
    double *wQuad_tmp, *xQuad_tmp;
    wQuad_tmp = new double[freq_points];
    xQuad_tmp = new double[freq_points];

    // -- Gauss-Leg quadrature -- //
    if(quad_type == 0)
    {
        gauleg(a,b,xQuad_tmp,wQuad_tmp,freq_points); // returns grid on [-1,1] from NR3
    }
    else if(quad_type == 1)
    {
        ReimannQuad(a,b,xQuad_tmp,wQuad_tmp,freq_points);
    }
    else
    {
        fprintf(stdout,"quadrature rule not coded\n");
        delete[] wQuad_tmp;
        delete[] xQuad_tmp;
        exit(1);
    }

    /* --- make weights of type gsl_complex_vector --- */
    gsl_complex zM1;
    for (int i = 0; i < freq_points; i++) 
    {
        GSL_SET_COMPLEX(&zM1,wQuad_tmp[i],0.0);
        gsl_vector_complex_set(wQuad_c,i,zM1);
    }

    for (int i = 0; i < freq_points; i++){
        gsl_vector_set(xQuad_c,i,xQuad_tmp[i]);
    }


    //std::cout << " frequency points \n";
    //OutputArray(freq_points,xQuad_tmp);
    //std::cout << "quad weights \n";
    //OutputArray(freq_points,wQuad_tmp);
    
    delete[] wQuad_tmp;
    delete[] xQuad_tmp;
}

void gsl_vector_sqrt(gsl_vector_complex *out, const gsl_vector_complex *in)
{

    gsl_complex z1;

    for(int i=0; i < in->size;i++)
    {
         z1 = gsl_complex_sqrt(gsl_vector_complex_get(in,i));
         gsl_vector_complex_set(out,i,z1);
    }

}

gsl_complex WeightedInner(const gsl_vector_complex *weights, const gsl_vector_complex *u, const gsl_vector_complex *v)
{
    // compute <u,v> = \sum_i (u_i^* v_i)*w_i ... weights w //

    gsl_complex ans;

    gsl_vector_complex *vc_tmp;
    vc_tmp = gsl_vector_complex_alloc(v->size);

    gsl_vector_complex_memcpy(vc_tmp,v);
    gsl_vector_complex_mul(vc_tmp,weights); // pointwise multiplcation, vc_tmp is updated
    gsl_blas_zdotc(u,vc_tmp,&ans); 

    //gsl_blas_zdotc(u,v,&ans);      // without weights this is MUCH faster
   
    gsl_vector_complex_free(vc_tmp);

    return ans;
}

double GetNorm_double(const gsl_vector_complex *u, const gsl_vector_complex *wQuad)
{
    gsl_complex nrmc = WeightedInner(wQuad,u,u);
    double nrm = gsl_complex_abs(nrmc);
    nrm = sqrt(nrm);
    //GSL_SET_COMPLEX(&nrmc,1.0/nrm,0.0);
    return nrm;
}

void NormalizeVector(gsl_vector_complex *u, const gsl_vector_complex *wQuad)
{
    gsl_complex nrmc;// = GetNorm_complex(u,wQuad);
    double nrm = GetNorm_double(u,wQuad);
    GSL_SET_COMPLEX(&nrmc,1.0/nrm,0.0);
    gsl_vector_complex_scale(u,nrmc);
}

// normalize row vectors using weight w //
void NormalizeTS(gsl_matrix_complex *A, const gsl_vector_complex *w)
{

    const int rows = A->size1;
    const int cols = A->size2;

    std::cout << "rows (NTS) = " << rows << std::endl;
    std::cout << "freq points (NTS) = " << cols << std::endl;

    gsl_complex sum_c, inv_norm_c;
    gsl_vector_complex *ts_el;
    double norm;

    ts_el = gsl_vector_complex_alloc(cols);

    for(int i=0; i < rows; i++)
    {

        gsl_matrix_complex_get_row(ts_el,A,i);
        NormalizeVector(ts_el,w);
        gsl_matrix_complex_set_row(A,i,ts_el);

    }

    gsl_vector_complex_free(ts_el);

}

void TF2_FullWaveform(gsl_vector_complex *wv, double *params, const gsl_vector *xQuad, double amp, double PN)
{
    // parameter list such that (m1(param),m2(param)) is a unique point in parameter space
    double TS_r = 0.0;
    double TS_i = 0.0;
    gsl_complex zM;

    for(int cols = 0; cols < xQuad->size; cols++)
    {
        TF2_Waveform(TS_r, TS_i, params, xQuad->data[cols], amp, PN);
        GSL_SET_COMPLEX(&zM, TS_r, TS_i);
        gsl_vector_complex_set(wv,cols,zM);
    }

}


void FillTrainingSet(gsl_matrix_complex *TS_gsl, const gsl_vector *xQuad, const gsl_vector_complex *wQuad, const TrainSet ts, const int rank, const bool whiten)
{
    // TODO: changes to this entire routine needed for different approximants //

    fprintf(stdout,"Populating training set on proc %i...\n",rank);

    // parameter list such that (m1(param),m2(param)) is a unique point in parameter space
    gsl_vector_complex *wv;
    double *params;
    double PN;
    int start_ind, end_ind, global_i, matrix_size;
    int gsl_status;
    gsl_vector *asd;
    gsl_complex z;
    double a;

    wv = gsl_vector_complex_alloc(xQuad->size);
    // if whitening, load the ASD
    if (whiten)
    {
        fprintf(stdout, "Loading ASD to whiten...\n");
        asd = gsl_vector_alloc(xQuad->size);
        FILE *asdf = fopen("ASD.txt", "r");
        gsl_status = gsl_vector_fscanf(asdf, asd);
        fclose(asdf);
        if( gsl_status == GSL_EFAILED )
        {
            fprintf(stderr, "Error reading ASD.txt\n");
            exit(1);
        }
    }


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
    

    if(strcmp(ts.model,"TaylorF2_PN3pt5") == 0)
    {
        fprintf(stdout,"Using the TaylorF2 spa approximant to PN=3.5\n");
        PN        = 3.5;
        params    = new double[4]; // (m1,m2,tc,phi_c)
        params[2] = 0.0;  // dummy variable (tc in waveform generation)
        params[3] = 0.0;  // dummy variable (phi_c in waveform generation)

        for(int i = 0; i < matrix_size; i++)
        {
            global_i = start_ind + i;
            params[0] = ts.m1[global_i] * mass_to_sec;
            params[1] = ts.m2[global_i] * mass_to_sec;
            TF2_FullWaveform(wv,params,xQuad,1.0,PN);
            // divide by the ASD to whiten
            if(whiten)
            {
                for(int jj = 0; jj < xQuad->size; jj++)
                {
                    z = gsl_vector_complex_get(wv, jj);
                    a = gsl_vector_get(asd, jj);
                    /*if(i < 10 && jj < 10)
                    {
                        fprintf(stdout, "i: %i j: %i\n  a: %e  z before: %e + i %e", i, jj, a, GSL_REAL(z), GSL_IMAG(z));
                    }*/
                    z = gsl_complex_div_real(z, a);
                    /*if(i < 10 && jj < 10)
                    {
                        fprintf(stdout, "  z after: %e + i %e\n", GSL_REAL(z), GSL_IMAG(z));
                    }*/
                    gsl_vector_complex_set(wv, jj, z);
                }
            }
            gsl_matrix_complex_set_row(TS_gsl,i,wv);
        }

    }
    else
    {
        std::cout << "Approximant not supported!" << std::endl;
        exit(1);
    }    

    std::cout << "Training set size is " << ts.ts_size << std::endl;

    /* -- Normalize training space here -- */
    fprintf(stdout,"Normalizing training set...\n");
    NormalizeTS(TS_gsl,wQuad);

    delete[] params;
    gsl_vector_complex_free(wv);
    if(whiten){
        gsl_vector_free(asd);}
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

void WriteGreedySelections(const int dim_RB,const int *selected_rows,const TrainSet ts)
{
    FILE *data;
    char filename[] = "GreedyPoints.txt";
    data = fopen(filename,"w");
    for(int i = 0; i < dim_RB ; i++)
    {
        fprintf(data,"%1.14f %1.14f\n",ts.m1[selected_rows[i]]/mass_to_sec,ts.m2[selected_rows[i]]/mass_to_sec);
    }
    fclose(data);
}

void WriteGreedyError(const int dim_RB, const double *app_err)
{

    FILE *data;
    char filename[] = "ApproxErrors.txt";
    data = fopen(filename,"w");
    for(int i = 0; i < dim_RB ; i++)
    {
        fprintf(data,"%1.14f\n",app_err[i]);
    }
    fclose(data);

}

void WriteGreedyRB(const gsl_matrix_complex *RB_space, const gsl_matrix_complex *R_matrix)
{
    // --- Output waveform here for matlab comparisons --- //
    FILE *rb_data, *r_data;
    char rb_filename[] = "Basis.txt";
    char r_filename[]  = "R.txt";

    //--- write R and RB to file ---//
    rb_data = fopen(rb_filename,"w");
    gsl_matrix_complex_fprintf(rb_data,RB_space, "%1.15e");
    fclose(rb_data);
    r_data = fopen(r_filename,"w");
    gsl_matrix_complex_fprintf(r_data,R_matrix, "%1.15e");
    fclose(r_data);

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

void GreedyMaster(const int size, const int max_RB, const int seed,const gsl_vector_complex *wQuad,const double tol, const TrainSet ts)
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
        IMGS(ru,ortho_basis,RB_space,wQuad,dim_RB);
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
    WriteGreedyRB(RB_space, R_matrix);
    WriteGreedySelections(dim_RB,greedy_points,ts);
    WriteGreedyError(dim_RB,greedy_err);


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

void Greedy(const int seed,const int max_RB, const gsl_matrix_complex *A,const gsl_vector_complex *wQuad,const double tol,const TrainSet ts)
{
// Input: 
//          A: gsl matrix of solutions (each row is a solutions, cols are quadrature samples)
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
    gsl_matrix_complex_set(R_matrix,0,0,tmpc);  // compute to mpi routine

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
        //IMGS(ru,ortho_basis,RB_space,wQuad,dim_RB);
        MGS(ru,ortho_basis,RB_space,wQuad,dim_RB);
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
    WriteGreedyRB(RB_space, R_matrix);
    WriteGreedySelections(dim_RB,greedy_points,ts);
    WriteGreedyError(dim_RB,greedy_err);


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

int main (int argc, char **argv) {

    // NOTE: use "new" for dynamic memory allocation, will be input in future version
    MPI::Init(argc, argv);


    // -- these information should be user input -- //
    // TODO: these can be input from a parameter file //
    const double a        = 40.0;                  // lower frequency value
    const double b        = 1024.0;                // upper frequency value
    const int freq_points = 1968;                  // total number of frequency points
    const int quad_type   = 1;                     // 0 = LGL, 1 = Reimman sum
    const int m_size      = 30;                    // parameter points in each m1,m2 direction
    const int ts_size     = 2734;                  // if reading ts from file, specify size
    const double m_low    = 1.0*mass_to_sec;       // lower mass value
    const double m_high   = 3.0*mass_to_sec;       // higher mass value
    bool load_from_file   = false;                 // load training points from file instead (file name used is below)
    const int seed        = 0;                     // greedy algorithm seed
    const double tol      = 1e-12;                 // greedy algorithm tolerance ( \| app \|^2)
    const char model_wv[] = "TaylorF2_PN3pt5";     // type of gravitational waveform model
    int max_RB            = 900;                   // estimated number of RB (reasonable upper bound)
    bool whiten           = false;                 // whether or not to whiten the waveforms with the ASD when calculating overlaps

    // -- declare variables --//
    int rank = 0;  // mpi information -- if running with mpi these will be reset
    int size = 1;  // mpi information -- if running with mpi these will be reset
    FILE *data1;
    gsl_matrix_complex *TS_gsl;
    gsl_vector_complex *wQuad;
    gsl_vector *xQuad;
    TrainSet ts;

    // Get number of procs this job is using (size) and unique rank of the processor this thread is running on //
    rank = MPI::COMM_WORLD.Get_rank();
    size = MPI::COMM_WORLD.Get_size();

    // Get the name of this processor (usually the hostname) //                        
    char name[MPI_MAX_PROCESSOR_NAME];
    int len;
    memset(name,0,MPI_MAX_PROCESSOR_NAME);
    MPI::Get_processor_name(name,len);
    memset(name+len,0,MPI_MAX_PROCESSOR_NAME-len);

    std::cout << "Number of tasks = " << size << " My rank = " << rank << " My name = " << name << "." << std::endl;


    // -- allocate memory --//
    wQuad = gsl_vector_complex_alloc(freq_points);
    xQuad = gsl_vector_alloc(freq_points);

    /* -- all procs should have a copy of the quadrature rule -- */
    // TODO: on nodes can this be shared? Whats the best way to impliment it? 
    MakeQuadratureRule(wQuad,xQuad,a,b,freq_points,quad_type);

    // -- build training set -- //
    // TODO: read in from file and populate ts //
    if(load_from_file)
    {
        BuildTS_from_file("test_bank.txt",ts_size,model_wv,ts);
    }
    else
    {
        BuildTS_tensor_product(m_size,m_low,m_high,model_wv,ts);
    }

    if(size == 1) // only 1 proc requested (serial mode)
    {
        TS_gsl = gsl_matrix_complex_alloc(ts.ts_size,freq_points); // GSL error handler will abort if too much requested
        FillTrainingSet(TS_gsl,xQuad,wQuad,ts,0,whiten);
        Greedy(seed,max_RB,TS_gsl,wQuad,tol,ts);
    }
    else
    {

        // -- split matrix TS_gsl among worker nodes. Assumes for-loop is "<" for this choice of myend -- //
        SplitTrainingSet(size,ts);

        if(rank != 0){
            TS_gsl = gsl_matrix_complex_alloc(ts.matrix_sub_size[rank-1],freq_points);
            FillTrainingSet(TS_gsl,xQuad,wQuad,ts,rank-1,whiten);
        }

        fprintf(stdout,"Finished distribution of training set\n");

        if(rank == 0){
            GreedyMaster(size,max_RB,seed,wQuad,tol,ts);
        }
        else{
            GreedyWorker(rank-1,max_RB,seed,wQuad,tol,TS_gsl,ts);
            gsl_matrix_complex_free(TS_gsl);
        }

    }

    if(rank == 0)
    {
        if(size == 1){
            WriteWaveform(xQuad->data,TS_gsl,0); // for comparison with other codes
            gsl_matrix_complex_free(TS_gsl);
        }

/*
        // --- output training set (assumes 1 proc) --- //
        char filename[] = "TS_Points.txt";
        data1 = fopen(filename,"w");
        for(int i = 0; i < ts.ts_size ; i++)
        {
            //fprintf(data1,"%.15lf %.15lf\n",ts.m1[i]/mass_to_sec,ts.m2[i]/mass_to_sec);
            fprintf(data1,"%.15le %.15le\n",ts.m1[i],ts.m2[i]);
        }
        fclose(data1);
*/

    }

    gsl_vector_complex_free(wQuad);
    gsl_vector_free(xQuad);
    free(ts.m1);
    free(ts.m2);

    if(ts.distributed){
        free(ts.mystart);
        free(ts.myend);
        free(ts.matrix_sub_size);
    }

    // Tell the MPI library to release all resources it is using
    MPI::Finalize();

}

