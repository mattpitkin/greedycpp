// AUTHOR :  Scott Field 
//           sfield@astro.cornell.edu
//
// DATE: Oct 1, 2014
//
// PURPOSE: various gsl helper routines

#ifndef gsl_helper_functions_h
#define gsl_helper_functions_h

#include <cmath>
#include <iostream>

#ifdef USE_NUMPY
//#include <complex.h>
#include "cnpy.h"
#include <complex>
#else
//#include <complex.h>
#include <cassert>
#include <iostream>
#include <complex>
#endif


// stick common gsl headers here -- this will reduce boilerplate. These gsl
// headers are used in every file which needed to include this hpp anyways.
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector_complex.h>

// Define OPTIMIZE_AVX with -D compiler flag // 
// compiling with icc and -O3 -xHOST flags should be faster

namespace mygsl {

#ifdef OPTIMIZE_AVX
// Returns the weighted complex scalar product u^H v working on the
// gsl_complex_vector data structures directly to avoid a copy,
// traversing the loop only once and enabling vectorization. 
gsl_complex WeightedInner(const gsl_vector_complex *weights,\
                          const gsl_vector_complex *u,\
                          const gsl_vector_complex *v);
#else
// --- compute <u,v> = \sum_i (u_i^* v_i)*w_i ... weights w --- //
// TODO: this routine can be made faster by fixing code!!
gsl_complex WeightedInner(const gsl_vector_complex *weights,
                          const gsl_vector_complex *u,
                          const gsl_vector_complex *v);
#endif

// --- compute <u,v> = \sum_i u_i^* v_i --- //
gsl_complex EuclideanInner(const gsl_vector_complex *u,
                           const gsl_vector_complex *v);

// --- compute <u,v> = \sum_i (u_i^* v_i)*w_i  --- //
// Will use faster routines if w_i is independent of i
gsl_complex InnerProduct(const gsl_vector_complex *weights,
                         const gsl_vector_complex *u,
                         const gsl_vector_complex *v,
                         const bool euclidean);

// tests whether a complex vector is constant
bool IsConstantVector(const gsl_vector_complex *u);

// --- returns WeightedInner(wQuad,u,u) as type double --- //
double GetNorm_double(const gsl_vector_complex *u,
                      const gsl_vector_complex *wQuad);

void NormalizeVector(gsl_vector_complex *u,
                     const gsl_vector_complex *wQuad);

// normalize row vectors of matrix A using weight w //
void NormalizeMatrixRows(gsl_matrix_complex *A, const gsl_vector_complex *w);

// Returns the index of the maximum value in the vector u.
// When there are several equal maximum
// elements then the lowest index is returned.
int gsl_vector_complex_max_index(const gsl_vector_complex *u);

void gsl_vector_complex_conj(gsl_vector_complex *u);

void gsl_vector_complex_parts(double *v_real,
                              double *v_imag,
                              const gsl_vector_complex * v_gsl);

void gsl_vector_complex_gsl_parts(gsl_vector *v_real,
                                  gsl_vector *v_imag,
                                  const gsl_vector_complex * v_gsl);

void make_gsl_vector_complex_parts(const double *v_real,
                                   const double *v_imag,
                                   gsl_vector_complex * v_gsl);

/* write real or imaginary part of m_gsl to file */
void gsl_matrix_complex_fprintf_part(const char *data_filename,
                                     const gsl_matrix_complex * m_gsl,
                                     const char *part);

/* write complex matrix m_gsl to file */
void gsl_matrix_complex_fprintf(const char *data_filename,
                                const gsl_matrix_complex * m_gsl);

// Lower triangular version of gsl_linalg_R_solve
int gsl_linalg_complex_L_solve(const gsl_matrix_complex* L,
                             const gsl_vector_complex* b,
                             gsl_vector_complex* x);

// sum first n rows of A's column i
// XXX: DUE TO MULTIPLE gets, THIS IS VERY SLOW
double SumColumn(const gsl_matrix_complex *A,
                           const int i,
                           const int n);

void gsl_vector_sqrt(gsl_vector_complex *out,
                     const gsl_vector_complex *in);

/* write complex matrix m_gsl to numpy file */
void gsl_matrix_complex_npy_save(const char *data_filename,
                                 const gsl_matrix_complex * m_gsl);

/* load complex matrix m_gsl from numpy file */
// m_gsl must be the correct shape
void gsl_matrix_complex_npy_load(const char *data_filename,
                                 gsl_matrix_complex * m_gsl);

// Modified GS routine //
void MGS(gsl_vector_complex *ru,
         gsl_vector_complex *ortho_basis,
         const gsl_matrix_complex *RB_space,
         const gsl_vector_complex *wQuad,
         const int dim_RB);

// Iterated modified GS routine //
void IMGS(gsl_vector_complex *ru,
          gsl_vector_complex *ortho_basis,
          const gsl_matrix_complex *RB_space,
          const gsl_vector_complex *wQuad,
          const int dim_RB);

}; // namespace mygsl

#endif /* gsl_helper_functions.h */
