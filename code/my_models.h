// AUTHOR :  Scott Field 
//           sfield@astro.cornell.edu
//
// DATE: Oct 1, 2014
//
// PURPOSE: model specific part of code. Edit as indicated below.

// *** THIS IS THE ONLY MODEL SPECIFIC PART OF THE CODE *** //


#ifndef my_models_h
#define my_models_h

#include "training_set.hpp"
#include "gsl_helper_functions.hpp"
#include <utmpx.h>

// model specific header files //
#include "../models/taylorf2/spa_waveforms.h"

#ifdef MODEL_LAL
#include "../models/lal/phenomp.h"
#include "../models/lal/SEOBNRv2_ROM.h"
#include "../models/lal/TaylorF2.h"
#endif

namespace mymodel {


// *** BEGIN MODEL SPECIFIC SECTION *** //
void EvaluateModel(gsl_vector_complex *model_eval,
                     const gsl_vector *xQuad,
                     const double *params,
                     const TrainingSetClass &ts)
{

  // New models go here...add to the list and loop over paramters //
  std::string model_tag(ts.model()); // TODO: ts.model should be std::string
  if(strcmp(ts.model(),"TaylorF2_PN3pt5") == 0) {
    TF2_FullWaveform(model_eval,params,xQuad,1.0,3.5); //amp=1.0,PN=3.5
  }
  #ifdef MODEL_LAL
  else if( strcmp(ts.model(),"PhenomP_plus")  == 0 ||
           strcmp(ts.model(),"PhenomP_cross") == 0 ||
           strcmp(ts.model(),"PhenomP_hphp")  == 0 ||
           strcmp(ts.model(),"PhenomP_hchc")  == 0 ||
           strcmp(ts.model(),"PhenomP_hphc")  == 0 ||
           strcmp(ts.model(),"PhenomP_all_parts") == 0) {
    PhenP_Waveform(model_eval, xQuad, params, model_tag);
  }
  else if (strcmp(ts.model(),"SEOBNRv2_ROM_SingleSpin") == 0)
    SEOBNRv2_ROM_SingleSpin_Waveform(model_eval, xQuad, params);
  else if (strcmp(ts.model(),"SEOBNRv2_ROM_DoubleSpin_HI") == 0)
    ROM_SEOBNRv2_DS_HI_FullWaveform(model_eval, xQuad, params);
  else if ( strcmp(ts.model(),"LackeyTidal2013_SEOBNRv2_ROM_HI_plus")  == 0 ||
            strcmp(ts.model(),"LackeyTidal2013_SEOBNRv2_ROM_HI_cross") == 0 ||
            strcmp(ts.model(),"LackeyTidal2013_SEOBNRv2_ROM_HI_hphp")  == 0 ||
            strcmp(ts.model(),"LackeyTidal2013_SEOBNRv2_ROM_HI_hchc")  == 0 ||
            strcmp(ts.model(),"LackeyTidal2013_SEOBNRv2_ROM_HI_hphc")  == 0 ||
            strcmp(ts.model(),"LackeyTidal2013_SEOBNRv2_ROM_HI_all_parts") == 0 ) {
    LackeyTidal2013_FullWaveform(model_eval, xQuad, params, model_tag);
  }
  else if (strcmp(ts.model(),"TaylorF2_LAL") == 0)
    TaylorF2_LAL_Waveform(model_eval, xQuad, params);
  #endif
  else {
    std::cerr << "my_models.h: Model not supported! Add model tag."<<std::endl;
    exit(1);
  }
}
// *** END MODEL SPECIFIC SECTION *** //


void FillTrainingSet(gsl_matrix_complex *TS_gsl,
                     const gsl_vector *xQuad,
                     const gsl_vector_complex *wQuad,
                     const TrainingSetClass &ts,
                     const int rank)
{

  if(rank <=40) {
    fprintf(stdout,"Populating training set on proc %i. Using model %s.\n",
            rank,ts.model());
  }

  int proc_ts_size;
  ts.LocalTrainingSetSize(proc_ts_size,rank);


  #ifdef USE_OPENMP // due to extra allocs, avoid this code if not using omp
  #pragma omp parallel
  {
    /*#pragma omp master
    {
      std::cout<<"threads (Fill TS matrix) = "<<omp_get_num_threads()<<std::endl;
    }*/

    //std::ostringstream os;
    //fprintf(stdout, "\nThread %i on cpu %i\n",omp_get_thread_num(),sched_getcpu());
    //std::cout<<os.str()<<std::flush;
    //std::cout<<num;

    double *params;
    gsl_vector_complex *model_eval;
    params     = new double[ts.param_dim()];
    model_eval = gsl_vector_complex_alloc(xQuad->size);
 
    #pragma omp for
    for(int i = 0; i < proc_ts_size; i++){
      ts.GetParameterValue(params,rank,i); //params [global_i][j] * (scale[j])
      EvaluateModel(model_eval,xQuad,params,ts);
      gsl_matrix_complex_set_row(TS_gsl,i,model_eval);
    }

    delete[] params;
    gsl_vector_complex_free(model_eval);

  }
  #else

  double *params;
  gsl_vector_complex *model_eval;
  params     = new double[ts.param_dim()];
  model_eval = gsl_vector_complex_alloc(xQuad->size);

  for(int i = 0; i < proc_ts_size; i++){
    ts.GetParameterValue(params,rank,i); //params [global_i][j] * (scale[j])
    EvaluateModel(model_eval,xQuad,params,ts);
    gsl_matrix_complex_set_row(TS_gsl,i,model_eval);
  }

  delete[] params;
  gsl_vector_complex_free(model_eval);

  #endif

  // -- Normalize training space here -- //
  //fprintf(stdout,"Normalizing training set...\n");
  mygsl::NormalizeMatrixRows(TS_gsl,wQuad); //TODO: not parallelized for OpenMP

}

};

#endif /* my_models.h */
