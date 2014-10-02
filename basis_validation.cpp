// AUTHOR :  Scott Field 
//           sfield@astro.cornell.edu
//
// DATE: Oct 1, 2014
//
// PURPOSE: test basis accuracy from random model samples

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <cmath>

#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_complex.h>
#include <gsl/gsl_block_complex_float.h>

#include "parameters.hpp"
#include "gsl_helper_functions.h"
#include "quadratures.h"
#include "training_set.hpp"
#include "my_models.h"

int main (int argc, char **argv) {

  //----- Checking the number of Variables passed to the Executable -----//
  if (argc != 2) {
    std::cerr << "Argument: 1. location of a cfg configuration/parameter file (ends in .cfg)" << std::endl;
    std::cerr << "Argument: 2. directory containing basis and quadrature information" << std::endl;
    //exit(0);
  }
  std::cout << "parameter file is: " << argv[1] << std::endl;
  std::cout << "basis folder is: " << argv[2] << std::endl;

  //--- Read input file argv[1]. If there is an error, report and exit.
  //--- Parameters class contains relevant information about parameters 
  Parameters *params_from_file = new Parameters(argv);
  std::cout << *params_from_file;

  gsl_matrix_complex *RB_space, *model_evaluations;
  gsl_vector_complex *wQuad, *model_eval, *r_tmp;
  gsl_vector *xQuad;
  double * errors;
  char err_filename[100];
  FILE *err_data;

  model_eval = gsl_vector_complex_alloc(params_from_file->quad_points());  
  r_tmp      = gsl_vector_complex_alloc(params_from_file->max_RB());
  RB_space   = gsl_matrix_complex_alloc(params_from_file->max_RB(),params_from_file->quad_points());

  // read in basis from gsl binary file //
  char rb_filename[100];
  FILE *pBASIS;
  strcpy(rb_filename,argv[2]);
  strcat(rb_filename,"/Basis.bin");
  std::cout << "reading basis from file " << rb_filename << " ..." << std::endl;
  pBASIS = fopen(rb_filename,"rb");
  gsl_matrix_complex_fread(pBASIS,RB_space);
  fclose(pBASIS);
  //char testwrite[100];
  //strcpy(testwrite,"RB_test");
  //mygsl::gsl_matrix_complex_fprintf(testwrite,RB_space);

  // reconstruct quadrature rule used in greedy //
  SetupQuadratureRule(&wQuad,&xQuad,params_from_file);

  // TODO: overload so random points can be filled here //
  TrainingSetClass *random_samples = new TrainingSetClass(params_from_file,1);

  model_evaluations = gsl_matrix_complex_alloc(random_samples->ts_size(),xQuad->size);
  errors            = new double[random_samples->ts_size()];

  mymodel::FillTrainingSet(model_evaluations,xQuad,wQuad,*random_samples,0);

  // TODO: use openMP for this part //
  // error reported will be \sqrt(h - Ph) //
  for(int ii = 0; ii < random_samples->ts_size(); ii++) {
    gsl_matrix_complex_get_row(model_eval,model_evaluations,ii);
    mygsl::MGS(r_tmp,model_eval,RB_space,wQuad,params_from_file->max_RB()-1);
    errors[ii] = GSL_REAL(gsl_vector_complex_get(r_tmp,params_from_file->max_RB()-1));
  }


  strcpy(err_filename,argv[2]);
  strcat(err_filename,"/VerErrors.txt");
  err_data = fopen(err_filename,"w");
  for(int i = 0; i < random_samples->ts_size() ; i++) {
    fprintf(err_data,"%1.15e\n",errors[i]);
  }
  fclose(err_data);


  gsl_matrix_complex_free(RB_space);
  gsl_matrix_complex_free(model_evaluations);
  gsl_vector_complex_free(r_tmp);
  gsl_vector_complex_free(model_eval);
  gsl_vector_complex_free(wQuad);
  gsl_vector_free(xQuad);

  delete [] errors;

  delete params_from_file;
  params_from_file = NULL;
}

