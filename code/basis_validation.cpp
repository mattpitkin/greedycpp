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

#ifdef USE_OPENMP
#include <omp.h>
#endif

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
#include "utils.h"

#include <assert.h>

int main (int argc, char **argv) {


  // -- record parameters with errors above tolerance to separate file -- //
  double err_tol = 1e-4;  // error recorded as \sqrt(h - Ph) 

  //----- Checking the number of Variables passed to the Executable -----//
  if (argc != 4) {
    std::cerr << 
         "Argument: 1. location of a cfg and parameter file" << std::endl;
    std::cerr << 
         "Argument: 2. directory (ending with /) with basis+quadrature info\n";
    std::cerr << "Argument: 3. file containing random samples" << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "parameter file is: " << argv[1] << std::endl;
  std::cout << "basis folder is: " << argv[2] << std::endl;
  std::cout << "random file is: " << argv[3] << std::endl;

  std::string random_sample_file = std::string(argv[3]);

  //--- Read input file argv[1]. If there is an error, report and exit.
  //--- Parameters class contains relevant information about parameters 
  Parameters *params_from_file = new Parameters(argv,true);
  std::cout << *params_from_file;

  gsl_matrix_complex *RB_space, *model_evaluations;
  gsl_vector_complex *wQuad;
  gsl_vector *xQuad;
  double * errors;
  char err_filename[200];
  char bad_param_filename[200];
  char shell_command[200];
  FILE *err_data, *bad_param;
  clock_t start, end;
  char rb_size_filename[200];
  char quad_size_filename[200];

  // -- Deduce quadrature and basis sizes -- //
  strcpy(rb_size_filename,argv[2]);
  strcat(rb_size_filename,"/ApproxErrors.txt");
  strcpy(quad_size_filename,argv[2]);
  strcat(quad_size_filename,"/quad_rule.txt");
  const int rb_size   = fcount_pts(rb_size_filename);
  const int quad_size = fcount_pts(quad_size_filename);

  RB_space = gsl_matrix_complex_alloc(rb_size,quad_size);

  // read in basis from gsl binary file //
  char rb_filename[200];
  FILE *pBASIS;
  strcpy(rb_filename,argv[2]);
  strcat(rb_filename,"/Basis.bin");
  std::cout << "reading basis from file " 
            << rb_filename << " ..." << std::endl;
  pBASIS = fopen(rb_filename,"rb");
  gsl_matrix_complex_fread(pBASIS,RB_space);
  fclose(pBASIS);

  // reconstruct quadrature rule used in greedy //
  SetupQuadratureRule(&wQuad,&xQuad,params_from_file);

  assert(quad_size==xQuad->size && quad_size==params_from_file->quad_points());

  // this is useful sanity check (looks at TS waveform errors) //
  //TrainingSetClass *random_samples = new TrainingSetClass(params_from_file,1);
  TrainingSetClass *random_samples = 
    new TrainingSetClass(params_from_file,random_sample_file);

  // Creating Run Directory //
  strcpy(shell_command, "mkdir -p -m700 ");
  strcat(shell_command, argv[2]);
  strcat(shell_command, "/validations/");
  system(shell_command);

  model_evaluations = 
   gsl_matrix_complex_alloc(random_samples->ts_size(),xQuad->size);
  errors            = new double[random_samples->ts_size()];

  mymodel::FillTrainingSet(model_evaluations,xQuad,wQuad,*random_samples,0);

  // error reported will be \sqrt(h - Ph) //
  start = clock();
  #ifdef USE_OPENMP
  double omp_start = omp_get_wtime();
  #endif

  #pragma omp parallel
  {
    // every variable declared here is thread private (thread-safe)
    gsl_vector_complex *model_eval, *r_tmp;
    model_eval = gsl_vector_complex_alloc(params_from_file->quad_points());  
    r_tmp      = gsl_vector_complex_alloc(rb_size);
 
    #pragma omp for
    for(int ii = 0; ii < random_samples->ts_size(); ii++) {
      gsl_matrix_complex_get_row(model_eval,model_evaluations,ii);
      mygsl::MGS(r_tmp,model_eval,RB_space,wQuad,rb_size-1);
      errors[ii] = GSL_REAL(gsl_vector_complex_get(r_tmp,rb_size-1));
      fprintf(stdout,"Random point index %i with error %1.3e\n",ii,errors[ii]);
    }

    gsl_vector_complex_free(r_tmp);
    gsl_vector_complex_free(model_eval);

  }
  
  end = clock();
  double alg_time = ((double) (end - start)/CLOCKS_PER_SEC);

  #ifdef USE_OPENMP
  double omp_end  = omp_get_wtime();
  double omp_time = omp_end - omp_start;
  #else
  double omp_time = alg_time;
  #endif

  // -- record results -- //
  fprintf(stdout,"validation took %f cpu secs and %f wall secs \n",
          alg_time,omp_time);
  strcpy(err_filename,argv[2]);
  strcat(err_filename,"validations/");
  strcpy(bad_param_filename,err_filename);
  strcat(bad_param_filename,"points_");
  strcat(err_filename,random_sample_file.substr(
    random_sample_file.find_last_of("\\/")+1,100).c_str());
  strcat(bad_param_filename,random_sample_file.substr(
    random_sample_file.find_last_of("\\/")+1,100).c_str());

  err_data  = fopen(err_filename,"w");
  bad_param = fopen(bad_param_filename,"w");
  for(int i = 0; i < random_samples->ts_size() ; i++) {
    random_samples->fprintf_ith(err_data,i);
    fprintf(err_data," %1.15e\n",errors[i]);
    if(errors[i]>err_tol) {
      random_samples->fprintf_ith(bad_param,i);
      fprintf(bad_param,"\n");
    }
  }
  fclose(bad_param);
  fclose(err_data);


  gsl_matrix_complex_free(RB_space);
  gsl_matrix_complex_free(model_evaluations);
  gsl_vector_complex_free(wQuad);
  gsl_vector_free(xQuad);

  delete [] errors;

  delete params_from_file;
  params_from_file = NULL;
}

