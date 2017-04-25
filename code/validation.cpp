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

int main (int argc, char **argv) {


  // -- record parameters with errors above tolerance to separate file -- //
  // TODO: read in tol from *.cfg file (should be fudge_factor*tol)
  //double err_tol = 1e-4;  // error recorded as \sqrt(h - Ph) --> this is projection (NOT eim)
  double err_tol = 1e-6;
    
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
  double *errors, *errors_eim, *model_norm;
  clock_t start, end;


  // If basis is orthogonal (but not normalized) carry out normalization //
  // NOTE: validation studies assume the basis satisfies <ei,ej> = \delta_{ij}
  //std::cout << "Normalizing the basis..." << std::endl;
  //mygsl::NormalizeMatrixRows(RB_space,wQuad);

  // -- this is useful sanity check (looks at TS waveform errors) -- //
  //TrainingSetClass *random_samples=new TrainingSetClass(&params_from_file,1);

  // -- use random samples file -- //
  TrainingSetClass *random_samples =
    new TrainingSetClass(&params_from_file,random_sample_file);

  // Use this if filling up matrix upfront (faster for smaller studies) 
  /*model_evaluations = 
   gsl_matrix_complex_alloc(random_samples->ts_size(),params_from_file.quad_points());
  mymodel::FillTrainingSet(model_evaluations,&xQuad,&wQuad,*random_samples,0);*/

  errors     = new double[random_samples->ts_size()];
  errors_eim = new double[random_samples->ts_size()];
  if(save_norm) {
    model_norm = new double[random_samples->ts_size()];
  }


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

      // Use this if model_evaluations matrix has been filled (see above) //
      //gsl_matrix_complex_get_row(model_eval,model_evaluations,ii);

      // Use this if model evaluations are done on-the-fly //
      //start1 = clock();
      random_samples->GetParameterValue(params,0,ii);
      mymodel::EvaluateModel(model_eval,&xQuad,params,*random_samples);

      // NOTE: not the best way to do this: recomputing the norm, nrm2 NOT squared.
      if(save_norm) {
        nrm = mygsl::GetNorm_double(model_eval,&wQuad);
        model_norm[ii] = nrm;
      }
      mygsl::NormalizeVector(model_eval,&wQuad);
      const double nrm2 = mygsl::GetNorm_double(model_eval,&wQuad); // if model = 0, norm = 0, else 1
      //fprintf(stdout,"nrm=%1.16e\n",nrm);

      /*end1 = clock();
      alg_time1 = ((double) (end1 - start1)/CLOCKS_PER_SEC);
      fprintf(stdout,"evaluating the model took %f seconds\n",alg_time1);*/

      // -- Compute empirical interpolation errors -- //
      // projection error validation modifies model_eval; 
      // so EIM errors must be computed first
      //start1 = clock();
      gsl_vector_complex *eim_eval =
        eim->eim_full_vector(model_eval,&RB_space,data->rb_size());
      gsl_vector_complex_sub(eim_eval,model_eval);

      // TODO: unaware of norm (compare with projection error below)
      errors_eim[ii] = mygsl::GetNorm_double(eim_eval,&wQuad);

      gsl_vector_complex_free(eim_eval);
      /*end1 = clock();
      alg_time1 = ((double) (end1 - start1)/CLOCKS_PER_SEC);
      fprintf(stdout,"eim took %f seconds\n",alg_time1);*/

      // -- Compute errors by projecting onto the basis -- //
      // NOTE: slow and fast way agree to about 1 sig fig
      //       MGS takes a more "stable" projection
      //start1 = clock();
      // slow way (older)
      //mygsl::MGS(r_tmp,model_eval,&RB_space,&wQuad,data->rb_size());
      //errors[ii] = GSL_REAL(gsl_vector_complex_get(r_tmp,data->rb_size()));
      // fast way (newer -- 6/15)
      gsl_vector_complex_mul(model_eval,&wQuad);
      mygsl::gsl_vector_complex_conj(model_eval);
      gsl_blas_zgemv(CblasNoTrans,GSL_COMPLEX_ONE,&RB_space,
                 model_eval,GSL_COMPLEX_ZERO,r_tmp);
      // to inspect projection coefficients
      /*for(int jj=0;jj<data->rb_size()+1;++jj) {
        std::cout << "size of r_tmp[jj] with jj = " << jj << " is = "
                  << gsl_complex_abs(gsl_vector_complex_get(r_tmp,jj)) 
                  << std::endl;
      }*/
      const double r_tmp_nrm = gsl_blas_dznrm2(r_tmp);
      double err_sqrd = nrm2 - r_tmp_nrm*r_tmp_nrm;
      if(err_sqrd < 0.0) // floating point error can trigger this
        errors[ii] = 1.0e-8;
      else
        errors[ii] = sqrt(err_sqrd);

      /*end1 = clock();
      alg_time1 = ((double) (end1 - start1)/CLOCKS_PER_SEC);
      fprintf(stdout,"MGS took %f seconds\n",alg_time1);*/
      //fprintf(stdout,"Random point index %i with error %1.3e\n",
      //  ii,errors[ii]);

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
  strcpy(model_norm_filename,err_filename);
  strcat(model_norm_filename,"model_norm_");

  strcat(err_filename,random_sample_file.substr(
    random_sample_file.find_last_of("\\/")+1,100).c_str());
  strcat(bad_param_filename,random_sample_file.substr(
    random_sample_file.find_last_of("\\/")+1,100).c_str());
  strcat(model_norm_filename,random_sample_file.substr(
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

  FILE *model_nrm = fopen(model_norm_filename,"w");
  if (bad_param==NULL) {
    std::cerr << "could not open model norm file\n";
    std::cerr << "with file name = " << model_norm_filename << std::endl;
    exit(1);
  }

  for(int i = 0; i < random_samples->ts_size() ; i++) {
    random_samples->fprintf_ith(err_data,i);
    fprintf(err_data," %1.15e",errors_eim[i]);
    fprintf(err_data," %1.15e\n",errors[i]);
    if(save_norm) {
      fprintf(model_nrm," %1.15e\n",model_norm[i]);
    }
    if(errors[i]>err_tol) {
      random_samples->fprintf_ith(bad_param,i);
      fprintf(bad_param,"\n");
    }
  }

  fclose(bad_param);
  fclose(err_data);
  fclose(model_nrm);

  if(save_norm) {
    delete [] model_norm;
  }
  delete [] errors;
  delete [] errors_eim;
  delete data;
  data = NULL;
  delete eim;
  eim = NULL;

  return 0;

}

