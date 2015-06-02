#ifndef quadratures_hpp
#define quadratures_hpp


class Parameters;

#include <gsl/gsl_vector_complex.h>
#include <gsl/gsl_vector_complex_double.h>
#include <gsl/gsl_matrix_complex_double.h>

/* --- fill array with linear spacing --- */
void Linspace(const int &n,const double &a,const double &b,double *SomeArray);

void ReimannQuad(const double a, const double b, double *xQuad, double * wQuad,
                 const int quad_points);

void MakeQuadratureRule(gsl_vector_complex *wQuad_c, 
                        gsl_vector *xQuad_c, 
                        const double a, 
                        const double b, 
                        const int quad_points,
                        const int quad_type);

void MakeWeightedInnerProduct(gsl_vector_complex *wQuad, FILE *weightf);

void SetupQuadratureRule(gsl_vector_complex **wQuad, 
                         gsl_vector **xQuad,
                         const Parameters *pParams);

#endif // quadratures.hpp //
