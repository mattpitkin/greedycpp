// AUTHOR :  Scott Field 
//           sfield@astro.cornell.edu
//
// DATE: Oct 1, 2014
//
// PURPOSE: various gsl helper routines


#ifndef gsl_helper_functions_h
#define gsl_helper_functions_h
#include <assert.h>

#ifdef USE_NUMPY
#include <complex.h>
#include "cnpy.h"
#include <complex>
#endif


// Define OPTIMIZE_AVX with -D compiler flag // 
// compiling with icc and -O3 -xHOST flags should be faster

namespace mygsl {

#ifdef OPTIMIZE_AVX
// Returns the weighted complex scalar product u^H v working on the
// gsl_complex_vector data structures directly to avoid a copy,
// traversing the loop only once and enabling vectorization. 
gsl_complex WeightedInner(const gsl_vector_complex *weights,\
                          const gsl_vector_complex *u,\
                          const gsl_vector_complex *v)
{
  gsl_complex ans;

  const size_t N = weights->size;

  assert(u->size == N && v->size == N);
    
  double ar = 0.0;
  double ai = 0.0;

  for (size_t i = 0; i < N; i++)
  {
    double weightsr = weights->data[2*i];
    double weightsi = weights->data[2*i+1];
    double ur = u->data[2*i];
    double ui = -u->data[2*i+1];
    double vr = v->data[2*i];
    double vi = v->data[2*i+1];
    double tmpr = weightsr*vr - weightsi*vi;
    double tmpi = weightsr*vi + weightsi*vr;
    ar += ur*tmpr - ui*tmpi;
    ai += ur*tmpi + ui*tmpr;
  }
  GSL_SET_COMPLEX(&ans, ar, ai);

  return ans;
}

#else

// --- compute <u,v> = \sum_i (u_i^* v_i)*w_i ... weights w --- //
// TODO: this routine can be made faster by fixing code!!
gsl_complex WeightedInner(const gsl_vector_complex *weights,\
                          const gsl_vector_complex *u,\
                          const gsl_vector_complex *v)
{
  gsl_complex ans;
  gsl_vector_complex *v_tmp;

  v_tmp = gsl_vector_complex_alloc(v->size);

  gsl_vector_complex_memcpy(v_tmp,v);
  gsl_vector_complex_mul(v_tmp,weights); // pointwise multiplcation
  gsl_blas_zdotc(u,v_tmp,&ans); 

  gsl_vector_complex_free(v_tmp);

  return ans;
}

#endif

// --- compute <u,v> = \sum_i u_i^* v_i --- //
gsl_complex EuclideanInner(const gsl_vector_complex *u,\
                           const gsl_vector_complex *v)
{
  gsl_complex ans;
  gsl_blas_zdotc(u,v,&ans);

  return ans;
}

// --- returns WeightedInner(wQuad,u,u) as type double --- //
double GetNorm_double(const gsl_vector_complex *u,\
                      const gsl_vector_complex *wQuad)
{
  gsl_complex nrmc = WeightedInner(wQuad,u,u);
  double nrm = gsl_complex_abs(nrmc);
  nrm = sqrt(nrm);
  return nrm;
}

void NormalizeVector(gsl_vector_complex *u,\
                     const gsl_vector_complex *wQuad)
{
  gsl_complex nrmc;
  double nrm = GetNorm_double(u,wQuad);
  GSL_SET_COMPLEX(&nrmc,1.0/nrm,0.0);
  gsl_vector_complex_scale(u,nrmc);
}

// normalize row vectors of matrix A using weight w //
void NormalizeMatrixRows(gsl_matrix_complex *A, const gsl_vector_complex *w)
{

  const int rows = A->size1;
  const int cols = A->size2;

  //std::cout << "rows of A (Training space) = " << rows << std::endl;
  //std::cout << "cols of A (Training space) = " << cols << std::endl;

  gsl_complex sum_c, inv_norm_c;
  gsl_vector_complex *row_vector;
  double norm;

  row_vector = gsl_vector_complex_alloc(cols);

  for(int i = 0; i < rows; i++) {
    gsl_matrix_complex_get_row(row_vector,A,i);
    NormalizeVector(row_vector,w);
    gsl_matrix_complex_set_row(A,i,row_vector);
  }

  gsl_vector_complex_free(row_vector);

}

void gsl_vector_complex_parts(double *v_real,\
                              double *v_imag,\
                              const gsl_vector_complex * v_gsl)
{
  for(int i = 0; i < v_gsl->size; i++) {
    v_real[i] = GSL_REAL(gsl_vector_complex_get(v_gsl,i));
    v_imag[i] = GSL_IMAG(gsl_vector_complex_get(v_gsl,i));
  }
}

void gsl_vector_complex_gsl_parts(gsl_vector *v_real,\
                                  gsl_vector *v_imag,\
                                  const gsl_vector_complex * v_gsl)
{

  double tmp_r, tmp_i;

  for(int i = 0; i < v_gsl->size; i++)
  {
    tmp_r = GSL_REAL(gsl_vector_complex_get(v_gsl,i));
    tmp_i = GSL_IMAG(gsl_vector_complex_get(v_gsl,i));

    gsl_vector_set(v_real, i,tmp_r);
    gsl_vector_set(v_imag, i,tmp_i);
  }
}

void make_gsl_vector_complex_parts(const double *v_real,\
                                   const double *v_imag,\
                                   gsl_vector_complex * v_gsl)
{
  gsl_complex tmpc;

  for(int i = 0; i < v_gsl->size; i++){
    GSL_SET_COMPLEX(&tmpc,v_real[i],v_imag[i]);
    gsl_vector_complex_set(v_gsl,i,tmpc);
  }
}

/* write real or imaginary part of m_gsl to file */
void gsl_matrix_complex_fprintf_part(const char *data_filename,\
                                     const gsl_matrix_complex * m_gsl,\
                                     const char *part)
{
  // NOTE: output so that A = R * Q is desired QR decomposition //

  const int cols = m_gsl->size2;
  const int rows = m_gsl->size1;

  FILE *pFILE;
  pFILE = fopen(data_filename,"w");

  gsl_vector_complex *v_gsl;
  double *vec_real, *vec_imag;
  v_gsl         = gsl_vector_complex_alloc(cols);
  vec_real      = (double *)malloc(cols*sizeof(double));
  vec_imag      = (double *)malloc(cols*sizeof(double));

  for(int ii = 0; ii < rows; ii++)
  {  
    gsl_matrix_complex_get_row(v_gsl,m_gsl, ii);
    gsl_vector_complex_parts(vec_real,vec_imag,v_gsl);

    for(int jj = 0; jj < cols; jj++)
    {  
      if(strcmp(part,"real") == 0){
        fprintf(pFILE, "%1.12e\t",vec_real[jj]);
      }
      else if(strcmp(part,"imag") == 0){
        fprintf(pFILE, "%1.12e\t",vec_imag[jj]);
      }
      else{
        fprintf(stderr,"part unknown");
        exit(1);
      }
    }

   fprintf(pFILE,"\n");

 }

  fclose(pFILE);

  gsl_vector_complex_free(v_gsl);
  free(vec_real);
  free(vec_imag);

}

/* write complex matrix m_gsl to file */
void gsl_matrix_complex_fprintf(const char *data_filename,\
                                     const gsl_matrix_complex * m_gsl)
{

  char data_filename_r[100];
  char data_filename_i[100];

  strcpy(data_filename_r,data_filename);
  strcat(data_filename_r,"_real.txt");
  strcpy(data_filename_i,data_filename);
  strcat(data_filename_i,"_imag.txt");

  gsl_matrix_complex_fprintf_part(data_filename_r,m_gsl,"real"); 
  gsl_matrix_complex_fprintf_part(data_filename_i,m_gsl,"imag"); 
}


// sum first n rows of A's column i
double SumColumn(const gsl_matrix_complex *A,
                           const int i,
                           const int n)
{
  double tmp = 0;
  gsl_complex tmpc;

  for(int j = 0; j < n; j++) {
    tmpc = gsl_matrix_complex_get(A,j,i);
    tmp += tmpc.dat[0]*tmpc.dat[0]+tmpc.dat[1]*tmpc.dat[1];
  }

  return tmp;
}


void gsl_vector_sqrt(gsl_vector_complex *out,\
                     const gsl_vector_complex *in)
{  
  gsl_complex z1;

  for(int i=0; i < in->size;i++){   
     z1 = gsl_complex_sqrt(gsl_vector_complex_get(in,i));
     gsl_vector_complex_set(out,i,z1);
  }
}

/* write complex matrix m_gsl to numpy file */
void gsl_matrix_complex_npy_save(const char *data_filename,\
                                 const gsl_matrix_complex * m_gsl)
{

  #ifdef USE_NUMPY
  const int cols = m_gsl->size2; // "quad_size"
  const int rows = m_gsl->size1; // "dim_rb"

  // TODO: this is going to double the memory footprint! //
  std::complex<double>* data = new std::complex<double>[cols*rows];
  gsl_vector_complex *row_vec;
  row_vec = gsl_vector_complex_alloc(cols);
  for(int ii=0;ii<rows;ii++) {
    gsl_matrix_complex_get_row(row_vec,m_gsl,ii);
    for(int jj=0;jj<cols;jj++) {
      data[cols*ii+jj] =  std::complex<double>(row_vec->data[2*jj]   + I*row_vec->data[2*jj+1]);
    }
  }

  const unsigned int shape[] = {rows,cols};
  cnpy::npy_save(data_filename,data,shape,2,"w");

  gsl_vector_complex_free(row_vec);
  delete [] data;
  #else
  std::cerr << "\nCode not compiled with cnpy\n";
  std::terminate();
  #endif

}

/* load complex matrix m_gsl from numpy file */
// m_gsl must be the correct shape
void gsl_matrix_complex_npy_load(const char *data_filename,\
                                 gsl_matrix_complex * m_gsl)
{

  #ifdef USE_NUMPY
  const int cols = m_gsl->size2; // "quad_size"
  const int rows = m_gsl->size1; // "dim_rb"

  // TODO: Large memory footprint here
  cnpy::NpyArray arr = cnpy::npy_load(data_filename);
  std::complex<double>* loaded_data =
    reinterpret_cast<std::complex<double>*>(arr.data);
  gsl_complex tmp_gsl;
  for(int ii=0;ii<rows;ii++) {
    for(int jj=0;jj<cols;jj++) {
      std::complex<double> tmp = loaded_data[cols*ii+jj];
      GSL_SET_COMPLEX(&tmp_gsl,tmp.real(),tmp.imag());
      gsl_matrix_complex_set(m_gsl, ii, jj, tmp_gsl);
    }
  }
  delete [] loaded_data;
  //arr.destruct();
  #else
  std::cerr << "\nCode not compiled with cnpy\n";
  std::terminate();
  #endif

}

void MGS(gsl_vector_complex *ru,\
         gsl_vector_complex *ortho_basis,\
         const gsl_matrix_complex *RB_space,\
         const gsl_vector_complex *wQuad,\
         const int dim_RB)
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

  // if not done, R matrix fills up below diagonal with .5 instead of 0
  gsl_vector_complex_set_zero(ru); 

  for(int i = 0; i < dim_RB; i++)
  {
    gsl_matrix_complex_get_row(basis,RB_space,i);

    /* --- ortho_basis = ortho_basis - L2_proj*basis; --- */
    L2_proj = WeightedInner(wQuad,basis,ortho_basis);
    gsl_vector_complex_set(ru,i,L2_proj);
    gsl_vector_complex_scale(basis,L2_proj); // basis <- basis*L2_proj
    gsl_vector_complex_sub(ortho_basis,basis); // ortho <- ortho - basis

  }

  double nrm = GetNorm_double(ortho_basis,wQuad);
  gsl_complex nrmc;
  GSL_SET_COMPLEX(&nrmc,nrm,0.0);
  gsl_vector_complex_set(ru,dim_RB,nrmc);

  NormalizeVector(ortho_basis,wQuad);

  gsl_vector_complex_free(basis);

}


void IMGS(gsl_vector_complex *ru,\
          gsl_vector_complex *ortho_basis,\
          const gsl_matrix_complex *RB_space,\
          const gsl_vector_complex *wQuad,\
          const int dim_RB)
{
/*  Iterated modified GS routine. 

   Input:  RB_space:    an existing orthonormal set of basis vectors
           ortho_basis: basis we shall orthogonalize agains RB_space
           wQuad:       quadrature weights for inner products
           dim_RB:      elements in RB_space (ortho_basis is dim_RB+1 element)
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

  // if not done, R matrix fills up below diagonal with .5 instead of 0
  gsl_vector_complex_set_zero(ru);  

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


    if( nrm_current/nrm_prev <= ortho_condition ) {
      nrm_prev = nrm_current;
      iter = iter + 1;
      gsl_vector_complex_memcpy(e,ortho_basis);
    }
    else{
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

}; // namespace mygsl

#endif /* gsl_helper_functions.h */
