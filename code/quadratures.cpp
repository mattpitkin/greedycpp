
#include "parameters.hpp"
#include "quadratures.hpp"
#include "gsl_helper_functions.hpp"

/* --- fill array with linear spacing --- */
void Linspace(const int &n, const double &a, const double &b, double *SomeArray)
{
  double factor = (b-a)/(double)(n-1);
  for(int i=0;i<n;i++) {
      SomeArray[i] = a + (double)i*factor;
  }
}

void ReimannQuad(const double a, const double b, double *xQuad, double * wQuad,
                 const int quad_points)
{
  Linspace(quad_points, a, b, xQuad);

  for(int i = 0; i < quad_points; i++) {
      wQuad[i] = (b-a)/( (double) (quad_points-1) );
  }

}

void MakeQuadratureRule(gsl_vector_complex *wQuad_c, 
                        gsl_vector *xQuad_c, 
                        const double a, 
                        const double b, 
                        const int quad_points,
                        const int quad_type)
{
  double *wQuad_tmp, *xQuad_tmp;
  wQuad_tmp = new double[quad_points];
  xQuad_tmp = new double[quad_points];

  // -- Quadrature rule for inner product between rows -- //
  if(quad_type == 0) {
      // TODO: add in gsl quad routine
      //gauleg(a,b,xQuad_tmp,wQuad_tmp,quad_points); // returns grid on [-1,1]
      fprintf(stderr,"gaussian quadrature rule removed\n");
  }
  else if(quad_type == 1) {
    ReimannQuad(a,b,xQuad_tmp,wQuad_tmp,quad_points);
  } // quad_type == 2 is handled before this routine is called
  else {
    fprintf(stderr,"quadrature rule not coded\n");
    delete[] wQuad_tmp;
    delete[] xQuad_tmp;
    exit(1);
  }

  /* --- Make quad weights and points of type gsl_complex_vector --- */
  gsl_complex zM1;
  for (int i = 0; i < quad_points; i++) {
    GSL_SET_COMPLEX(&zM1,wQuad_tmp[i],0.0);
    gsl_vector_complex_set(wQuad_c,i,zM1);
  }

  if(quad_type != 2) {
    for (int i = 0; i < quad_points; i++){
      gsl_vector_set(xQuad_c,i,xQuad_tmp[i]);
     }
  }

  delete[] wQuad_tmp;
  delete[] xQuad_tmp;
}

void MakeWeightedInnerProduct(gsl_vector_complex *wQuad, FILE *weightf)
{
/* If the inner product includes a weight W(x), modify the quadrature weights.
   routine called if 'weighted = true' from configuration (input) file. */

  // TODO: error check that input weight and wQuad are of equal length //

  int gsl_status;
  gsl_complex z;
  double a;
  gsl_vector *weight;

  // load the weight //
  fprintf(stdout, "Loading ASD (weight) ...\n");
  weight = gsl_vector_alloc(wQuad->size);
  gsl_status = gsl_vector_fscanf(weightf, weight);
  if( gsl_status == GSL_EFAILED ) {
    fprintf(stderr, "Error reading ASD file\n");
    exit(1);
  }

  // divide by the ASD //
  for(int jj = 0; jj < wQuad->size; jj++) {
    z = gsl_vector_complex_get(wQuad, jj);
    a = gsl_vector_get(weight, jj);

    // TODO: should this be squared for ASD?
    z = gsl_complex_div_real(z, a);

    gsl_vector_complex_set(wQuad, jj, z);
  }

    gsl_vector_free(weight);

}

void SetupQuadratureRule(gsl_vector_complex **wQuad, 
                         gsl_vector **xQuad,
                         const Parameters *pParams)
{
    // wQuad and xQuad are pointers to pointers, free memory in main

  const int quad_type = pParams->quad_type();
  const bool weighted_inner = pParams->weighted_inner();
  const int quad_points = pParams->quad_points();
  const double x_min = pParams->x_min();
  const double x_max = pParams->x_max();
  const char *quad_nodes_file = strdup(pParams->quad_nodes_file().c_str());
  const char *weight_file_name = strdup(pParams->quad_weight_file().c_str());
  const char *num_weight_file = strdup(pParams->num_weight_file().c_str());

  int gsl_status;

  gsl_vector_complex *wQuad_tmp;
  gsl_vector *xQuad_tmp, *wQuad_tmp1;

  // -- allocate memory --//
  wQuad_tmp  = gsl_vector_complex_alloc(quad_points);
  xQuad_tmp  = gsl_vector_alloc(quad_points);
  wQuad_tmp1 = gsl_vector_alloc(quad_points);
  gsl_complex zM1;

  // -- load quadrature nodes from file -- //
  if(quad_type == 2) {
    FILE *fvecf = fopen(quad_nodes_file, "r");
    if (fvecf==NULL) {
      fprintf(stderr,"Could not open quadrature nodes file \"%s\".\n", quad_nodes_file);
      exit(1);
    }
    gsl_status = gsl_vector_fscanf(fvecf, xQuad_tmp);
    fclose(fvecf);
    if( gsl_status == GSL_EFAILED ){
      fprintf(stderr, "Error reading qaudrature file %s\n",quad_nodes_file);
      exit(1);
    }

    FILE *fweightf = fopen(num_weight_file, "r");
    if (fweightf==NULL) {
      fprintf(stderr,"Could not open quadrature weights file \"%s\".\n", weight_file_name);
      exit(1);
    }
    gsl_status = gsl_vector_fscanf(fweightf, wQuad_tmp1);
    fclose(fweightf);
    if( gsl_status == GSL_EFAILED ){
      fprintf(stderr, "Error: numerical weights from %s\n", num_weight_file);
      exit(1);
    }

    for (int i = 0; i < quad_points; i++) {
      GSL_SET_COMPLEX(&zM1,gsl_vector_get(wQuad_tmp1,i),0.0);
      gsl_vector_complex_set(wQuad_tmp,i,zM1);
    }

  }
  else{
    /* -- all procs should have a copy of the quadrature rule -- */
    // TODO: on nodes can this be shared? Whats the best way to impliment it? 
    MakeQuadratureRule(wQuad_tmp,xQuad_tmp,x_min,x_max,quad_points,quad_type);
  }

  // Role inner product weight into wQuad //
  if(weighted_inner){
    FILE *weightf = fopen(weight_file_name, "r");
    MakeWeightedInnerProduct(wQuad_tmp, weightf);
    fclose(weightf);
  }

  *wQuad = wQuad_tmp;
  *xQuad = xQuad_tmp;

  gsl_vector_free(wQuad_tmp1);

}


