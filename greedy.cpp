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


/*-- constants --*/
double c = 299792458.0;                             // speed of light
double c2 = c*c;
double mass_to_sec = 1.98892e30*6.67428E-11/(c2*c); // Conversion factor to change Solar Mass units to seconds. assumes G = 6.67428E-11 and SolarMass = 1.98892e30
double DL = 3.08568025E22;                          // Luminosity distance in meters, set to one megaparsec


/*-- Training Set data structure -- */
typedef struct
tagTrainSet
{
  double *m1, *m2;      // each element of the training set is (m1[i],m2[i])
  int ts_size;          // number of training set elements
  bool distributed;     // true if training set will be split over worker procs/nodes
  int *matrix_sub_size; // ts.myend[rank] - ts.mystart[rank] --> number of rows on each proc
  int *mystart, *myend; // maps global row index of A into local worker index (when distributed = true)
  char model[100];      // name of gravitational waveform model (ex: TaylorF2_PN3pt5)
}
TrainSet;


/* --- fill array with linear spacing --- */
void Linspace(const int &n, const double &a, const double &b, double *SomeArray)
{
    double factor = (b-a)/(double)(n-1);
    for(int i=0;i<n;i++)
    {
        SomeArray[i] = a + (double)i*factor;
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

void BuildTS(const int &m_size, const double &m_low, const double &m_high, const char *model_name, TrainSet &ts)
{

    double *m1_tmp, *m2_tmp, *mass_list;
    double m_i, m_j;
    int counter = 0;

    m1_tmp = (double *)malloc(m_size*m_size*sizeof(double));
    m2_tmp = (double *)malloc(m_size*m_size*sizeof(double));
    mass_list = (double *)malloc(m_size*sizeof(double));

    if(m1_tmp==NULL || m2_tmp==NULL || mass_list==NULL)
    {
        std::cout << "Failed to allocate memory in BuildTS" << std::endl;
        free(m1_tmp); free(m2_tmp); free(mass_list);
        exit(1);
    }

    Linspace(m_size, m_low, m_high, mass_list); // m_size equidistant points from m_low to m_high

    // mass_list x mass_list tensor product -- (m1_temp[i],m2_temp[i]) ith training set element
    for(int i = 0; i < m_size; i++)
    {
        m_i = mass_list[i];

        for(int j = 0; j < m_size; j++)
        {
            m_j = mass_list[j];
            m1_tmp[counter] = m_i;
            m2_tmp[counter] = m_j;
            counter = counter + 1;
        }
    }

    /* --- fill training set data structure --- */
    ts.ts_size = m_size*m_size;
    ts.m2 = m2_tmp;
    ts.m1 = m1_tmp;
    ts.distributed = false; // by default. set to true if SplitTrainingSet is called
    strcpy(ts.model,model_name);

    std::cout << "Using waveform model: " << ts.model << std::endl;

    // NOTE: DONT free m1_tmp, m2_tmp here... do in main through ts

    free(mass_list);

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

void NormalizeVector(gsl_vector_complex *u, const gsl_vector_complex *wQuad)
{
    gsl_complex nrmc = WeightedInner(wQuad,u,u);
    double nrm = gsl_complex_abs(nrmc);
    nrm = sqrt(nrm);
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


void FillTrainingSet(gsl_matrix_complex *TS_gsl, const gsl_vector *xQuad, const gsl_vector_complex *wQuad, const TrainSet ts, const int rank)
{
    // TODO: changes to this entire routine needed for different approximants //

    fprintf(stdout,"Populating training set on proc %i...\n",rank);

    // parameter list such that (m1(param),m2(param)) is a unique point in parameter space
    gsl_vector_complex *wv;
    double *params;
    double PN;
    int start_ind, end_ind, global_i, matrix_size;

    wv = gsl_vector_complex_alloc(xQuad->size);

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

        //for(int rows = start_ind; rows < end_ind; rows++)
        //{

        for(int i = 0; i < matrix_size; i++)
        {
            global_i = start_ind + i;

            //std::cout << "global i " << global_i << std::endl;
            //params[0] = ts.m1[rows];
            //params[1] = ts.m2[rows];
            params[0] = ts.m1[global_i];
            params[1] = ts.m2[global_i];
            TF2_FullWaveform(wv,params,xQuad,1.0,PN);
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
}

void MGS(gsl_vector_complex *ortho_basis,const gsl_matrix_complex *RB_space,const gsl_vector_complex *wQuad, const int dim_RB)
{

    int quad_num = RB_space->size2;
    gsl_complex L2_proj, tmp;
    gsl_vector_complex *basis;


    basis = gsl_vector_complex_alloc(quad_num);


    for(int i = 0; i < dim_RB; i++)
    {
        gsl_matrix_complex_get_row(basis,RB_space,i);

        /* --- ortho_basis = ortho_basis - L2_proj*basis; --- */
        L2_proj = WeightedInner(wQuad,basis,ortho_basis);
        gsl_vector_complex_scale(basis,L2_proj); // basis <- basis*L2_proj
        gsl_vector_complex_sub(ortho_basis,basis); // ortho_basis <- ortho_basis - basis

    }

    NormalizeVector(ortho_basis,wQuad);

    gsl_vector_complex_free(basis);

}

void SplitTrainingSet(const int size, TrainSet &ts)
{
    int *mystart_tmp, *myend_tmp, *matrix_sub_size_tmp;

    int numprocs_worker = size - 1; // master (proc 0) will do no work
    int rows = ts.ts_size;
    int proc_id;


    mystart_tmp         = (int *)malloc((size-2)*sizeof(int));
    myend_tmp           = (int *)malloc((size-2)*sizeof(int));
    matrix_sub_size_tmp = (int *)malloc((size-2)*sizeof(int));

    for(int i = 2; i <= size; i++) 
    {
        proc_id = i-2;
        mystart_tmp[i-2] = (rows / numprocs_worker) * proc_id;
        if (rows % numprocs_worker > proc_id)
        {
            mystart_tmp[i-2] += proc_id;
            myend_tmp[i-2] = mystart_tmp[i-2] + (rows / numprocs_worker) + 1;
        }
        else
        {
            mystart_tmp[i-2] += rows % numprocs_worker;
            myend_tmp[i-2] = mystart_tmp[i-2] + (rows / numprocs_worker);
        }
        //fprintf(stdout,"Proc %i starting at row %i and ending at row %i",i-1,mystart_tmp[i-2],myend_tmp[i-2]);
    }

    ts.distributed = true;

    for(int i = 0; i < size - 1; i++){
        matrix_sub_size_tmp[i] = myend_tmp[i] - mystart_tmp[i];
    }

    // NOTE: DONT free m1_tmp, m2_tmp here, do this in main
    ts.mystart         = mystart_tmp;
    ts.myend           = myend_tmp;
    ts.matrix_sub_size = matrix_sub_size_tmp;

}

void GreedyWorker(const int rank, const gsl_vector_complex *wQuad,const double tol, const gsl_matrix_complex *A, const TrainSet ts)
{
/* worker routine for computing computationally intensive part of greedy */
    int est_RB = 312;          // estimated number of RB (reasonable upper bound)
    int dim_RB, worst_global, worst_local, pass_rb;
    double worst_app = 1.0;
    double cols = wQuad->size;
    double rb_inc_r, rb_inc_i, tmp;
    double *errors;
    gsl_complex tmpc;
    gsl_vector_complex *last_rb, *ts_el, *ortho_basis;
    gsl_matrix_complex *project_coeff;

    last_rb       = gsl_vector_complex_alloc(cols);
    ts_el         = gsl_vector_complex_alloc(cols);
    project_coeff = gsl_matrix_complex_alloc(est_RB,ts.matrix_sub_size[rank]);
    errors        = (double *)malloc(ts.matrix_sub_size[rank]*sizeof(double));
    ortho_basis   = gsl_vector_complex_alloc(cols);

    fprintf(stdout,"I'm worker %i and I was given %i matrix elements from %i to %i\n",rank,ts.matrix_sub_size[rank],ts.mystart[rank],ts.myend[rank]-1);


    int tmp_i;
    while(worst_app > tol)
    {
//        std:: cout << "proc " << rank << " is here 0" << std::endl;
        MPI_Barrier(MPI_COMM_WORLD); 
        MPI_Recv(&dim_RB, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);


        /* -- receive new rb -- */
        // TODO:this should be broadcasted for speed 
        for(int i = 0; i < cols; i++){
            MPI_Recv(&rb_inc_r, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            MPI_Recv(&rb_inc_i, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            GSL_SET_COMPLEX(&tmpc,rb_inc_r,rb_inc_i);
            gsl_vector_complex_set(last_rb,i,tmpc);
        }

//        std:: cout << "proc " << rank << " is here 1" << std::endl;

        // Compute overlaps of pieces of A with rb_new
        for(int i = 0; i < ts.matrix_sub_size[rank]; i++)
        {
            //tmp_i = ts.mystart[rank] + i; // TODO: fix this
            tmp_i = i;
            //fprintf(stdout,"worker %i global matrix index %i\n",rank,tmp_i);
            gsl_matrix_complex_get_row(ts_el,A,tmp_i);
            tmpc = WeightedInner(wQuad,last_rb,ts_el);
            gsl_matrix_complex_set(project_coeff,dim_RB-1,i,tmpc);

            tmp = 0;
            for(int j = 0; j < dim_RB; j++)
            {
                tmpc = gsl_matrix_complex_get(project_coeff,j,i);
                tmp = tmp + gsl_complex_abs(tmpc)*gsl_complex_abs(tmpc);
            }
            errors[i] = 1.0 - tmp;
        }

//        std:: cout << "proc " << rank << " is here 2" << std::endl;


        /* -- find worst error here -- */
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

//        std:: cout << "proc " << rank << " is here 3" << std::endl;


        /* -- pass worst error and index to master --*/
        MPI_Send(&worst_global, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&tmp, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

        /* -- pass row which is globally worst to master as new basis -- */
        //MPI_Recv(&pass_rb, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        // TODO: what if errors are equal?? which proc returns the necessary element
        // TODO: only pass necessary basis
        //gsl_matrix_complex_get_row(ortho_basis,A,worst_global);
        gsl_matrix_complex_get_row(ortho_basis,A,worst_local);
        for(int j = 0; j< cols; j++){
            MPI_Send(&GSL_REAL(gsl_vector_complex_get(ortho_basis,j)),1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&GSL_IMAG(gsl_vector_complex_get(ortho_basis,j)),1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }


//        std:: cout << "proc " << rank << " is here 4" << std::endl;

        MPI_Barrier(MPI_COMM_WORLD); // TODO:this should be broadcasted for speed 
        MPI_Recv(&worst_app, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }

    gsl_vector_complex_free(last_rb);
    gsl_vector_complex_free(ts_el);
    gsl_vector_complex_free(ortho_basis);
    gsl_matrix_complex_free(project_coeff);
    free(errors);
}

void GreedyMaster(const int size, const gsl_vector_complex *seed,const gsl_vector_complex *wQuad,const double tol, const TrainSet ts,int **sel_rows,double **app_err,int &dim_RB)
{
/* Input: 
          A: gsl matrix of solutions (each row is a solutions, cols are quadrature samples)
          seed: first greedy pick (the row itself)
          tol: approximation tolerance
          size: number of procs. if 1 serial mode assumed

   Output:
          sel_rows: row index defining reduced basis. sel_rows[0] = seed
          dim_RB: number of greedy_points
*/

    fprintf(stdout,"Starting greedy algorithm...\n");

    const int rows = ts.ts_size;  // number of rows to approximate
    const int cols = wQuad->size; // samples (for quadrature)
    int est_RB = 312;             // estimated number of RB (reasonable upper bound)
    int *greedy_points;           // selectted greedy points (row selection)
    double *greedy_err;           // approximate error
    clock_t start, end;           // for algorithm timing experiments
    double alg_time;              // for algorithm timing experiments
    double tmp,worst_err;         // errors in greedy sweep
    int worst_app, worst_worker, worst_rank, tmpi;             // worst error stored
    int *worst_workers;
    double *worst_errs;
    gsl_complex tmpc;          // worst error temp
    double rb_inc_r, rb_inc_i;

    /* --- Output waveform here for matlab comparisons --- */
    FILE *data;
    char filename[] = "Basis.txt";

    gsl_vector_complex *last_rb, *ortho_basis;  // for holding training space elements
    gsl_matrix_complex *RB_space;

    /* --- this memory should be freed here --- */
    last_rb       = gsl_vector_complex_alloc(cols);
    ortho_basis   = gsl_vector_complex_alloc(cols);
    worst_workers = (int *)malloc((size-1)*sizeof(int));
    worst_errs    = (double *)malloc((size-1)*sizeof(double));

    /* --- this memory should be freed in main --- */
    greedy_points = (int *)malloc(est_RB*sizeof(int));
    greedy_err    = (double *)malloc(est_RB*sizeof(double));
    RB_space      = gsl_matrix_complex_alloc(est_RB,cols); 

    /* --- initialize algorithm with seed --- */
    gsl_matrix_complex_set_row(RB_space,0,seed);
    greedy_points[0] = 0; // TODO: this should be set by input
    dim_RB           = 1;
    greedy_err[0]    = 1.0;


    /* --- Continue approximation until tolerance satisfied --- */
    start = clock();
    while(greedy_err[dim_RB-1] > tol)
    {
        gsl_matrix_complex_get_row(last_rb,RB_space,dim_RB-1); // get last computed basis

        MPI_Barrier(MPI_COMM_WORLD); // TODO:this should be broadcasted for speed

        for(int i = 1; i <= size-1; i++){
            MPI_Send(&dim_RB,1, MPI_INT, i, 0, MPI_COMM_WORLD);
            for(int j = 0; j< cols; j++){
                MPI_Send(&GSL_REAL(gsl_vector_complex_get(last_rb,j)),1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
                MPI_Send(&GSL_IMAG(gsl_vector_complex_get(last_rb,j)),1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
            }
        }


        for(int i = 0; i < size - 1; i++){
            MPI_Recv(&worst_workers[i], 1, MPI_INT, i+1, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            MPI_Recv(&worst_errs[i], 1, MPI_DOUBLE, i+1, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            //fprintf(stdout,"RB %i suggested by worker %i with error %1.14e\n",worst_workers[i],i,worst_errs[i]);
        }

        worst_err = 0.0;
        for(int i = 0; i < size - 1; i++){
            if(worst_err < worst_errs[i])
            {
                worst_err  = worst_errs[i];
                worst_app  = worst_workers[i];
                worst_rank = i+1;
            }
        }


        // TODO: this can be written a lot better (e.g. only pass the necesary one)
        for(int i = 1; i <= size-1; i++){

            /* -- receive suggestd next rb from each worker -- */
            for(int j = 0; j < cols; j++){
                MPI_Recv(&rb_inc_r, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                MPI_Recv(&rb_inc_i, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);

                if(i==worst_rank){
                    GSL_SET_COMPLEX(&tmpc,rb_inc_r,rb_inc_i);
                    gsl_vector_complex_set(ortho_basis,j,tmpc);
                }
            }

        }

        MPI_Barrier(MPI_COMM_WORLD);

        for(int i = 1; i <= size-1; i++){
            MPI_Send(&worst_err, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }

        /* --- add worst approximated element to basis --- */
        greedy_points[dim_RB] = worst_app;
        greedy_err[dim_RB] = worst_err;

        /* --- add worst approximated solution to basis set --- */
        MGS(ortho_basis,RB_space,wQuad,dim_RB);
        gsl_matrix_complex_set_row(RB_space,dim_RB,ortho_basis);
        dim_RB = dim_RB + 1;

        fprintf(stdout,"RB dimension %i || Current row selection %i || Approximation error %1.14e\n",dim_RB,worst_app,worst_err);

    }
    end = clock();

    alg_time = ((double) (end - start)/CLOCKS_PER_SEC);
    fprintf(stdout,"Building approximation space took %f seconds\n",alg_time);
    dim_RB = dim_RB - 1;
    *sel_rows = greedy_points;
    *app_err = greedy_err;

    gsl_vector_complex_free(last_rb);
    gsl_vector_complex_free(ortho_basis);
    gsl_matrix_complex_free(RB_space);
    free(worst_workers);
    free(worst_errs);    
}

void Greedy(const int seed,const gsl_matrix_complex *A,const gsl_vector_complex *wQuad,const double tol, const TrainSet ts,int **sel_rows,double **app_err,int &dim_RB)
{
/* Input: 
          A: gsl matrix of solutions (each row is a solutions, cols are quadrature samples)
          seed: first greedy pick
          tol: approximation tolerance

   Output:
          sel_rows: row index defining reduced basis. sel_rows[0] = seed
          dim_RB: number of greedy_points
*/

    fprintf(stdout,"Starting greedy algorithm in serial mode...\n");

    const int rows = A->size1; // number of rows to approximate
    const int cols = A->size2; // samples (for quadrature)
    int est_RB = 312;          // estimated number of RB (reasonable upper bound)
    int *greedy_points;        // selectted greedy points (row selection)
    double *greedy_err;        // approximate error
    clock_t start, end;        // for algorithm timing experiments
    double alg_time;           // for algorithm timing experiments
    double tmp,worst_err;      // errors in greedy sweep
    int worst_app;             // worst error stored
    gsl_complex tmpc;          // worst error temp


    /* --- Output waveform here for matlab comparisons --- */
    FILE *data;
    char filename[] = "Basis.txt";

    gsl_vector_complex *ts_el, *last_rb, *ortho_basis;  // for holding training space elements
    gsl_matrix_complex *RB_space;
    double *errors;                       // approximation errors with RB space of dimension dim_RB
    gsl_matrix_complex *project_coeff; // h = coeff_i e_i is approximation we seek

    /* --- this memory should be freed here --- */
    ts_el         = gsl_vector_complex_alloc(cols);
    last_rb       = gsl_vector_complex_alloc(cols);
    ortho_basis   = gsl_vector_complex_alloc(cols);
    errors        = (double *)malloc(rows*sizeof(double));
    project_coeff = gsl_matrix_complex_alloc(est_RB,rows);

    /* --- this memory should be freed in main --- */
    greedy_points = (int *)malloc(est_RB*sizeof(int));
    greedy_err    = (double *)malloc(est_RB*sizeof(double));
    RB_space      = gsl_matrix_complex_alloc(est_RB,cols); 

    /* --- initialize algorithm with seed --- */
    gsl_matrix_complex_get_row(ts_el,A,seed);
    gsl_matrix_complex_set_row(RB_space,0,ts_el);
    greedy_points[0] = seed;
    dim_RB           = 1;
    greedy_err[0]    = 1.0;

    /* --- Continue approximation until tolerance satisfied --- */
    start = clock();
    while(greedy_err[dim_RB-1] > tol)
    {

        gsl_matrix_complex_get_row(last_rb,RB_space,dim_RB-1); // get last computed basis
        worst_err = 0.0;

        /* --- Loop over training set ---*/
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

        /* --- add worst approximated element to basis --- */
        greedy_points[dim_RB] = worst_app;
        greedy_err[dim_RB] = worst_err;

        /* --- add worst approximated solution to basis set --- */
        gsl_matrix_complex_get_row(ortho_basis,A,worst_app);
        MGS(ortho_basis,RB_space,wQuad,dim_RB);
        gsl_matrix_complex_set_row(RB_space,dim_RB,ortho_basis);
        dim_RB = dim_RB + 1;

        fprintf(stdout,"RB dimension %i || Current row selection %i || Approximation error %1.14e\n",dim_RB,worst_app,worst_err);

    }
    end = clock();

    alg_time = ((double) (end - start)/CLOCKS_PER_SEC);
    fprintf(stdout,"Building approximation space took %f seconds\n",alg_time);
    dim_RB = dim_RB - 1;
    *sel_rows = greedy_points;
    *app_err = greedy_err;

    gsl_vector_complex_free(ts_el);
    gsl_vector_complex_free(last_rb);
    gsl_vector_complex_free(ortho_basis);
    free(errors);
    gsl_matrix_complex_free(project_coeff);
    gsl_matrix_complex_free(RB_space);
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


int main (int argc, char **argv) {

    // NOTE: use "new" for dynamic memory allocation, will be input in future version

    // -- these information should be user input -- //
    // TODO: these can be input from a parameter file //
    const double a        = 40.0;                  // lower frequency value
    const double b        = 1024.0;                // upper frequency value
    const int freq_points = 1969;                  // total number of frequency points
    const int quad_type   = 1;                     // 0 = LGL, 1 = Reimman sum
    const int m_size      = 30;                    // parameter points in each m1,m2 direction
    const double m_low    = 1.0*mass_to_sec;       // lower mass value
    const double m_high   = 3.0*mass_to_sec;       // higher mass value
    const int seed        = 0;                     // greedy algorithm seed
    const double tol      = 1e-12;                 // greedy algorithm tolerance ( \| app \|^2)
    const char model_wv[] = "TaylorF2_PN3pt5";     // type of gravitational waveform model

    // -- declare variables --//
    double *app_err;
    int *selected_rows;
    int dim_RB;
    int rank, size; // mpi information 
    gsl_matrix_complex *TS_gsl;
    gsl_vector_complex *wQuad, *seed_el;
    gsl_vector *xQuad;
    TrainSet ts;

    // -- allocate memory --//
    wQuad = gsl_vector_complex_alloc(freq_points);
    xQuad = gsl_vector_alloc(freq_points);

    /* -- all procs should have a copy of the quadrature rule -- */
    // TODO: on single node jobs can this be shared? Whats the best way to impliment it? 
    MakeQuadratureRule(wQuad,xQuad,a,b,freq_points,quad_type);

    // -- build training set -- //
    // TODO: read in from file and populate ts //
    BuildTS(m_size,m_low,m_high,model_wv,ts);



    // Initialize the MPI library, get number of procs this job is using (size) and unique rank of the processor this thread is running on 
    MPI::Init(argc, argv);
    rank = MPI::COMM_WORLD.Get_rank();
    size = MPI::COMM_WORLD.Get_size();

    // Get the name of this processor (usually the hostname).  We call                                                      
    // memset to ensure the string is null-terminated.  Not all MPI                                                        
    // implementations null-terminate the processor name since the MPI                                                     
    // standard specifies that the name is *not* supposed to be returned                                                   
    // null-terminated.                                                                                                    
    char name[MPI_MAX_PROCESSOR_NAME];
    int len;
    memset(name,0,MPI_MAX_PROCESSOR_NAME);
    MPI::Get_processor_name(name,len);
    memset(name+len,0,MPI_MAX_PROCESSOR_NAME-len);

    std::cout << "Number of tasks = " << size << " My rank = " << rank << " My name = " << name << "." << std::endl;

    double *params; // TODO hack for seed, fix
   
    if(size == 1) // only 1 proc requested (serial mode)
    {
        TS_gsl = gsl_matrix_complex_alloc(ts.ts_size,freq_points); // GSL error handler will abort if too much requested
        std::cout << "Im here 1 " << std::endl;

        FillTrainingSet(TS_gsl,xQuad,wQuad,ts,0);
        std::cout << "Im here 2 " << std::endl;

        Greedy(seed,TS_gsl,wQuad,tol,ts,&selected_rows,&app_err,dim_RB);
        std::cout << "Im here 4 " << std::endl;


        gsl_matrix_complex_free(TS_gsl);
    }
    else
    {

        /* -- split matrix TS_gsl among worker nodes. Assumes for-loop is "<" for this choice of myend -- */
        SplitTrainingSet(size,ts);

        //TS_gsl = gsl_matrix_complex_alloc(ts.ts_size,freq_points); // GSL error handler will abort if too much requested
        //FillTrainingSet(TS_gsl,xQuad,wQuad,ts,rank-1);
        //FillTrainingSet(TS_gsl,xQuad,wQuad,ts,0);

        //TODO: seed should be implimented in a better way
        if(rank == 0)
        {

            params = new double[4]; // (m1,m2,tc,phi_c)
            params[2] = 0.0;  // dummy variable (tc in waveform generation)
            params[3] = 0.0;  // dummy variable (phi_c in waveform generation)
            params[0] = ts.m1[seed];
            params[1] = ts.m2[seed];

            seed_el = gsl_vector_complex_alloc(freq_points);
            TF2_FullWaveform(seed_el,params,xQuad,1.0,3.5);
            NormalizeVector(seed_el,wQuad);

            delete[] params;
        }
        else
        {
            TS_gsl = gsl_matrix_complex_alloc(ts.matrix_sub_size[rank-1],freq_points); // GSL error handler will abort if too much requested
            FillTrainingSet(TS_gsl,xQuad,wQuad,ts,rank-1);

            //TS_gsl = gsl_matrix_complex_alloc(ts.ts_size,freq_points); // GSL error handler will abort if too much requested
            //FillTrainingSet(TS_gsl,xQuad,wQuad,ts,0);

        }

        MPI_Barrier(MPI_COMM_WORLD);
        fprintf(stdout,"Finished distribution of training set\n");

        if(rank == 0){

            //seed_el = gsl_vector_complex_alloc(freq_points);
            //gsl_matrix_complex_get_row(seed_el,TS_gsl,seed);

            GreedyMaster(size,seed_el,wQuad,tol,ts,&selected_rows,&app_err,dim_RB);

            //gsl_matrix_complex_free(TS_gsl);

        }
        else{
            GreedyWorker(rank-1,wQuad,tol,TS_gsl,ts);
            gsl_matrix_complex_free(TS_gsl);
        }

    }


    if(rank == 0)
    {
        /* --- output quantities of interest to file --- */
        WriteGreedySelections(dim_RB,selected_rows,ts);
        WriteGreedyError(dim_RB,app_err);
        //WriteWaveform(xQuad->data,TS_gsl,0); // for comparison with other codes

        free(app_err);
        free(selected_rows);
    }

    // Tell the MPI library to release all resources it is using:
    MPI::Finalize();

    gsl_vector_complex_free(wQuad);
    gsl_vector_free(xQuad);
    free(ts.m1);
    free(ts.m2);

    if(ts.distributed){
        //std::cout << "This was a dist run rank " << rank << std::endl;
        free(ts.mystart);
        free(ts.myend);
        free(ts.matrix_sub_size);
    }

}

