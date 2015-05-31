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

#ifdef USE_NUMPY
#include <complex.h>
#include "cnpy.h"
#include <complex>
#endif

#ifdef USE_OPENMP
#include <omp.h>
#endif

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector_complex.h>
#include <gsl/gsl_block_complex_float.h>

#include "load_simulation_data.hpp"
#include "parameters.hpp"
#include "gsl_helper_functions.hpp"
#include "quadratures.hpp"
#include "training_set.hpp"
#include "my_models.h"

int main (int argc, char **argv) {


  // -- record parameters with errors above tolerance to separate file -- //
  double err_tol = 1e-4;  // error recorded as \sqrt(h - Ph) 


  //----- Checking the number of Variables passed to the Executable -----//
  if (argc != 5) {
    std::cerr << 
         "Argument: 1. location of a cfg and parameter file" << std::endl;
    std::cerr << 
         "Argument: 2. directory (ending with /) with basis+quadrature info\n";
    std::cerr << "Argument: 3. file containing random samples" << std::endl;
    std::cerr << "Argument: 4. basis format (npy or gsl)" << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "parameter file is: " << argv[1] << std::endl;
  std::cout << "basis folder is: " << argv[2] << std::endl;
  std::cout << "random file is: " << argv[3] << std::endl;
  std::cout << "basis format is: " << argv[4] << std::endl;

  std::string basis_path         = std::string(argv[2]);
  std::string random_sample_file = std::string(argv[3]);

  if(basis_path.at(basis_path.length()-1) != '/') {
    std::cerr << "directory with basis file must end with '/' " << std::endl;
    exit(1);
  }

  //--- Load greedy data from previous simulation -- //
  LoadData *data = new LoadData(argv);
  const Parameters& params_from_file = data->params_from_file();
  const gsl_matrix_complex RB_space = data->RB_space();
  const gsl_vector_complex wQuad = data->wQuad();
  const gsl_vector xQuad = data->xQuad();

  gsl_matrix_complex *model_evaluations;
  double * errors;
  char err_filename[200];
  char bad_param_filename[200];
  char shell_command[200];
  FILE *err_data, *bad_param;
  clock_t start, end;


  // If basis is orthogonal (but not normalized) carry out normalization //
  // NOTE: validation studies assume the basis satisfies <ei,ej> = \delta_{ij}
  //std::cout << "Normalizing the basis..." << std::endl;
  //mygsl::NormalizeMatrixRows(RB_space,wQuad);

  // this is useful sanity check (looks at TS waveform errors) //
  //TrainingSetClass *random_samples = new TrainingSetClass(params_from_file,1);
  TrainingSetClass *random_samples = 
    new TrainingSetClass(&params_from_file,random_sample_file);

  // Creating Run Directory //
  strcpy(shell_command, "mkdir -p -m700 ");
  strcat(shell_command, argv[2]);
  strcat(shell_command, "/validations/");
  system(shell_command);

  // Use this if filling up matrix upfront (faster for smaller studies) 
  /*model_evaluations = 
   gsl_matrix_complex_alloc(random_samples->ts_size(),xQuad->size);
  mymodel::FillTrainingSet(model_evaluations,xQuad,wQuad,*random_samples,0);*/

  errors            = new double[random_samples->ts_size()];


  // error reported will be \sqrt(h - Ph) //
  start = clock();
  #ifdef USE_OPENMP
  double omp_start = omp_get_wtime();
  #endif

  #pragma omp parallel
  {

    // every variable declared here is thread private (thread-safe)
    gsl_vector_complex *model_eval, *r_tmp;
    model_eval = gsl_vector_complex_alloc(params_from_file.quad_points());  
    r_tmp      = gsl_vector_complex_alloc(data->rb_size());
    double *params;
    params     = new double[random_samples->param_dim()]; //

    #pragma omp for
    for(int ii = 0; ii < random_samples->ts_size(); ii++) {

      // Use this if model_evaluations matrix has been filled (see above) //
      //gsl_matrix_complex_get_row(model_eval,model_evaluations,ii);

      // Use this if model evalutions are done on-the-fly //
      random_samples->GetParameterValue(params,0,ii);
      mymodel::EvaluateModel(model_eval,&xQuad,params,*random_samples);
      mygsl::NormalizeVector(model_eval,&wQuad);

      mygsl::MGS(r_tmp,model_eval,&RB_space,&wQuad,data->rb_size()-1);
      errors[ii] = GSL_REAL(gsl_vector_complex_get(r_tmp,data->rb_size()-1));
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
  if (err_data==NULL) {
    std::cerr << "Could not open error report file\n";
    std::cerr << "with file name = " << err_filename << std::endl;
    exit(1);
  }

  bad_param = fopen(bad_param_filename,"w");
  if (bad_param==NULL) {
    std::cerr << "could not open bad param file\n";
    std::cerr << "with file name = " << bad_param_filename << std::endl;
    exit(1);
  }

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

  delete [] errors;
  delete data;
  data = NULL;
}

