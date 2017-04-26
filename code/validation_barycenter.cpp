// AUTHOR :  Scott Field 
//           sfield@astro.cornell.edu
//
// DATE: Oct 1, 2014
//
// PURPOSE: test basis accuracy from random model samples

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>

#ifdef USE_NUMPY
//#include <complex.h>
#include "cnpy.h"
#include <complex>
#endif

#ifdef USE_OPENMP
#include <omp.h>
#endif

#include "gsl_helper_functions.hpp"
#include "load_simulation_data.hpp"
#include "parameters.hpp"
#include "training_set.hpp"
#include "my_models.h"
#include "eim.hpp"

#define MAXABS(a,b)  fabs(a) > fabs(b) ? fabs(a) : fabs(b)
#define TWOPI 6.283185307179586476925286766559005768

int main (int argc, char **argv) {
  // -- record parameters with errors above tolerance to separate file -- //
  double err_tol = 1e-7; // maximum allowed residual time (seconds)
    
  bool save_norm = false;

  //----- Checking the number of Variables passed to the Executable -----//
  if (argc < 5 || argc > 6) {
    std::cerr << 
         "Argument: 1. location of a validation configuration (cfg) file\n";
    std::cerr << 
         "Argument: 2. directory (ending with /) with basis+quadrature info\n";
    std::cerr << "Argument: 3. basis format (npy or gsl)" << std::endl;
    std::cerr << "Argument: 4. file containing random samples" << std::endl;
    std::cerr << "Argument: 5. (Optional, without /) Local output directory [default=validations]\n";
    return EXIT_FAILURE;
  }

  std::string loc_dir("validations");
  if (argc == 6)
    loc_dir.assign(argv[5]);

  std::cout << "parameter file is: "      << argv[1] << std::endl;
  std::cout << "basis folder is: "        << argv[2] << std::endl;
  std::cout << "basis format is: "        << argv[3] << std::endl;
  std::cout << "random file is: "         << argv[4] << std::endl;
  std::cout << "local output directory: " << loc_dir << std::endl;


  std::string basis_path         = std::string(argv[2]);
  std::string random_sample_file = std::string(argv[4]);

  if(basis_path.at(basis_path.length()-1) != '/') {
    std::cerr << "directory with basis file must end with '/' " << std::endl;
    exit(1);
  }

  if(loc_dir.at(loc_dir.length()-1) == '/') {
    std::cerr << "output directory cannot end with '/' " << std::endl;
    exit(1);
  }

  char err_filename[200];
  char bad_param_filename[200];
  char model_norm_filename[200];
  char shell_command[200];

  // Creating Run Directory //
  loc_dir.append("/");
  loc_dir.insert(0,"/");
  std::cout << "loc dir is = " << loc_dir << std::endl;
  strcpy(shell_command, "mkdir -p -m700 ");
  strcat(shell_command, argv[2]);
  strcat(shell_command, loc_dir.c_str());
  int ret = system(shell_command);
  if(ret == -1) {
    std::cerr << "Could not make a run directory" << std::endl;
    exit(1);
  }

  // Copy validation run cfg file to the output folder //
  std::string validation_cfg(argv[2]);
  validation_cfg.append(loc_dir).append("validation_run.cfg");
  std::ifstream src(argv[1],std::ios::binary);
  std::ofstream dst(validation_cfg.c_str(),std::ios::binary);
  dst << src.rdbuf();
  src.close();
  dst.close();
  std::cout << "parameter file is: " << argv[1] << std::endl;

  //--- Load greedy data from previous simulation -- //
  LoadData *data = new LoadData(argv,true);
  EIM *eim       = new EIM(data->rb_size(),data->quad_size(),
                           data->invV(),data->eim_indices());

  const Parameters& params_from_file = data->params_from_file();
  const gsl_matrix_complex RB_space = data->RB_space();
  const gsl_vector_complex wQuad = data->wQuad();
  const gsl_vector xQuad = data->xQuad();

  gsl_matrix_complex *model_evaluations;
  
  // calculate the maximum residual time between the validation points and empirically interpolated points
  // and calculate the phase mismatch, assuming a 1000Hz source, caused by any barycentering uncertainty
  // 1 - sum(cos(2*pi*1000*(Dt)))/N
  double *error_resmax, *error_phasemismatch;
  clock_t start, end;

  // -- use random samples file -- //
  TrainingSetClass *random_samples =
    new TrainingSetClass(&params_from_file,random_sample_file);

  error_resmax        = new double[random_samples->ts_size()];
  error_phasemismatch = new double[random_samples->ts_size()];

  // error reported will be \sqrt(h - Ph) //
  start = clock();
  clock_t start1, end1;
  
  #ifdef USE_OPENMP
  double omp_start = omp_get_wtime();
  #endif

  #pragma omp parallel
  {

    // every variable declared here is thread private (thread-safe)
    gsl_vector_complex *model_eval, *r_tmp;
    model_eval = gsl_vector_complex_alloc(params_from_file.quad_points());

    // r_tmp is rb_size+1... in MGS routine this entry is ||ortho||
    // use this r_tmp for older way (see below)
    //r_tmp      = gsl_vector_complex_alloc(data->rb_size()+1);
    //gsl_vector_complex_set_zero(r_tmp);
    r_tmp      = gsl_vector_complex_alloc(data->rb_size());

    double *params;
    params     = new double[random_samples->param_dim()];
    double alg_time1;

    // estimate the size of loop allocated to each thread
    #ifdef USE_OPENMP
    int size_per_thread =
      std::floor(random_samples->ts_size()/omp_get_num_threads());
    //int one_percent_finished = std::floor(size_per_thread/100);
    int one_percent_finished = std::ceil(size_per_thread/100);
    if (one_percent_finished == 0) { // needs to be > 0
      one_percent_finished = 1;
    }
    int percent_completed = 0;
    int thread_id = omp_get_thread_num();
    fprintf(stdout,"Thread %i, estimates %i of %i total. %i for one percent\n",
            thread_id,size_per_thread,random_samples->ts_size(),
            one_percent_finished);
    #endif
    double nrm;

    #pragma omp for
    for(int ii = 0; ii < random_samples->ts_size(); ii++) {
      // Use this if model evaluations are done on-the-fly //
      random_samples->GetParameterValue(params,0,ii);
      mymodel::EvaluateModel(model_eval,&xQuad,params,*random_samples);

      // -- Compute empirical interpolation errors (residuals) -- //
      gsl_vector_complex *eim_eval = eim->eim_full_vector(model_eval,&RB_space,data->rb_size());
      
      //if ( ii == 0 ){
        // print out model and EIM model
        //FILE *fp1 = NULL;
        //fp1 = fopen("/home/matthew/testing/redordbar/validation/eim_model.txt", "w");
        //if ( gsl_vector_complex_fprintf(fp1, eim_eval, "%.16le") ){
        //  fprintf(stderr, "Problem outputting EIM model\n");
        //}
        //fclose(fp1);
        //fp1 = fopen("/home/matthew/testing/redordbar/validation/full_model.txt", "w");
        //if ( gsl_vector_complex_fprintf(fp1, model_eval, "%.16le") ){
        //  fprintf(stderr, "Problem outputting full model\n");
        //}
        //fclose(fp1);
      //}
        
      gsl_vector_complex_sub(eim_eval,model_eval);

      // get maximum absolute residual value (model is only real so get real vector view)
      gsl_vector_view realres = gsl_vector_complex_real(eim_eval);
      double resmin = 0., resmax = 0.;
      gsl_vector_minmax(&realres.vector, &resmin, &resmax); // get the minimum and maximum values
      error_resmax[ii] = MAXABS(resmin, resmax); // determine the maximum absolute value

      // get the phase mismatch for a 1000 Hz signal (use trapzeium rule for integration)
      double phasematch = 0.;
      for (int jj=0; jj < data->quad_size()-1; jj++){ phasematch += 0.5*(cos(TWOPI*1000.*gsl_vector_get(&realres.vector, jj)) + cos(TWOPI*1000.*gsl_vector_get(&realres.vector, jj+1))); }
      fprintf(stderr, "phase match = %.16lf\n", phasematch);
      phasematch /= (data->quad_size()-1.);
      error_phasemismatch[ii] = 1. - phasematch;

      gsl_vector_complex_free(eim_eval);

      #ifdef USE_OPENMP
      if( ii % one_percent_finished == 0) {
        percent_completed +=1;
        fprintf(stdout,"Thread %i %i percent finished\n",thread_id,
                percent_completed);
      }
      #endif
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

  // set err_filename and bad_param_filename //
  strcpy(err_filename,argv[2]);
  strcat(err_filename,loc_dir.c_str());
  strcpy(bad_param_filename,err_filename);
  strcat(bad_param_filename,"bad_points_");

  strcat(err_filename,random_sample_file.substr(
    random_sample_file.find_last_of("\\/")+1,100).c_str());
  strcat(bad_param_filename,random_sample_file.substr(
    random_sample_file.find_last_of("\\/")+1,100).c_str());

  FILE *err_data  = fopen(err_filename,"w");
  if (err_data==NULL) {
    std::cerr << "Could not open error report file\n";
    std::cerr << "with file name = " << err_filename << std::endl;
    exit(1);
  }

  FILE *bad_param = fopen(bad_param_filename,"w");
  if (bad_param==NULL) {
    std::cerr << "could not open bad param file\n";
    std::cerr << "with file name = " << bad_param_filename << std::endl;
    exit(1);
  }

  for(int i = 0; i < random_samples->ts_size() ; i++) {
    random_samples->fprintf_ith(err_data,i);
    fprintf(err_data," %1.15e",error_phasemismatch[i]);
    fprintf(err_data," %1.15e\n",error_resmax[i]);
    if(error_resmax[i]>err_tol) {
      random_samples->fprintf_ith(bad_param,i);
      fprintf(bad_param,"\n");
    }
  }

  fclose(bad_param);
  fclose(err_data);

  delete [] error_resmax;
  delete [] error_phasemismatch;
  delete data;
  data = NULL;
  delete eim;
  eim = NULL;

  return 0;

}

