// AUTHOR :  Scott Field 
//           sfield@umd.edu
//
// DATE: Jan 20, 2013
//
// PURPOSE: mpi version of greedy (pivoted MGS) algorithm 


//#include "nr3.h"
#include "gauss_wgts.h"

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
  double *m1, *m2; // each element of the training set is (m1[i],m2[i])
  int ts_size;     // number of training set elements
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


void BuildTS(const int &m_size, const double &m_low, const double &m_high, TrainSet &ts)
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
    // NOTE: DONT free m1_tmp, m2_tmp here... do in main through ts

    free(mass_list);

}

void OutputArray(const int n, double *list)
{
    for(int i=0;i<n;i++){
        std::cout << list[i] << std::endl;
    }
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

//void Greedy(const int seed,const gsl_matrix_complex *A,const gls_vector_complex *wQuad,const double tol,int *greedy_points,int dim_RB)
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
    gsl_complex ip, L2_proj;
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

            L2_proj = WeightedInner(wQuad,last_rb,ts_el);

            gsl_matrix_complex_set(project_coeff,dim_RB-1,i,L2_proj);

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
    const double a = 40.0;                  // lower frequency value
    const double b = 1024.0;                // upper frequency value
    const int freq_points = 1969;           // total number of frequency points
    const int m_size = 30;                  // parameter points in each m1,m2 direction
    const double m_low = 1.0*mass_to_sec;   // lower mass value
    const double m_high = 3.0*mass_to_sec;  // higher mass value
    const int seed = 0;                     // greedy algorithm seed
    const double tol = 1e-12;               // greedy algorithm tolerance ( \| app \|^2)
    const double GWamp = 1.0;               // GW amplitude
    const double PN = 3.5;                  // waveform's PN order

    // -- declare variables --//
    double *xQuad, *wQuad, *params;
    double *app_err;
    int *selected_rows;
    int dim_RB;
    gsl_matrix_complex *TS_gsl;
    gsl_vector_complex *wQuad_c;
    TrainSet ts;

    // -- build training set -- //
    BuildTS(m_size,m_low,m_high,ts);

    // -- allocate memory --//
    wQuad = new double[freq_points];
    xQuad = new double[freq_points];
    params = new double[4]; // (m1,m2,tc,phi_c)
    TS_gsl = gsl_matrix_complex_alloc(ts.ts_size,freq_points); // GSL error handler will abort if too much requested
    wQuad_c = gsl_vector_complex_alloc(freq_points);

    // -- Gauss-Leg quadrature -- //
    //gauleg(a,b,xQuad,wQuad,freq_points); // returns grid on [-1,1] from NR3
    ReimannQuad(a,b,xQuad,wQuad,freq_points);

    /* --- make weights of type gsl_complex_vector --- */
    gsl_complex zM1;
    for (int i = 0; i < freq_points; i++) 
    {
        GSL_SET_COMPLEX(&zM1,wQuad[i],0.0);
        gsl_vector_complex_set(wQuad_c,i,zM1);
    }


    //std::cout << " frequency points \n";
    //OutputArray(freq_points,xQuad);
    //std::cout << "quad weights \n";
    //OutputArray(freq_points,wQuad);
    

    // -- populate training space -- //
    // parameter list such that (m1(param),m2(param)) is a unique point in parameter space
    double TS_r = 0.0;
    double TS_i = 0.0;
    gsl_complex zM;

    fprintf(stdout,"Populating training set...\n");

    for(int rows = 0; rows < ts.ts_size; rows++)
    {

        params[0] = ts.m1[rows];
        params[1] = ts.m2[rows];
        params[2] = 0.0;  // dummy variable (tc in waveform generation)
        params[3] = 0.0;  // dummy variable (phi_c in waveform generation)

        for(int cols = 0; cols < freq_points; cols++)
        {

            TF2_Waveform(TS_r, TS_i, params, xQuad[cols], GWamp, PN);
            GSL_SET_COMPLEX(&zM, TS_r, TS_i);
            gsl_matrix_complex_set(TS_gsl,rows,cols,zM);

            // TO DO: Adding PSD should go here //
        }

    }

    std::cout << "Training set size is " << TS_gsl->size1 << std::endl;


    /* -- Normalize training space here -- */
    fprintf(stdout,"Normalizing training set...\n");
    NormalizeTS(TS_gsl,wQuad_c);

    /* -- greedy algorithm here, same in/out as MGSCP -- */
    fprintf(stdout,"Starting greedy algorithm...\n");
    Greedy(seed,TS_gsl,wQuad_c,tol,ts,&selected_rows,&app_err,dim_RB);

    /* --- output quantities of interest to file --- */
    WriteGreedySelections(dim_RB,selected_rows,ts);
    WriteGreedyError(dim_RB,app_err);
    WriteWaveform(xQuad,TS_gsl,0); // for comparison with other codes


    /* --- free from memory --- */
    delete[] xQuad;
    delete[] wQuad;
    delete[] params;

    gsl_matrix_complex_free(TS_gsl);
    gsl_vector_complex_free(wQuad_c);

    free(app_err);
    free(selected_rows);

    free(ts.m1);
    free(ts.m2);
}

