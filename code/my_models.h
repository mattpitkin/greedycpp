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

// model specific header files //
#include "../models/taylorf2/spa_waveforms.h"

#ifdef MODEL_LAL
#include "../models/lal/phenomp.h"
#include "../models/lal/SEOBNRv2_ROM.h"
#endif

namespace mymodel {


// *** BEGIN MODEL SPECIFIC SECTION *** //
void EvaluateModel(gsl_vector_complex *model_eval,
                     const gsl_vector *xQuad,
                     const double *params,
                     const TrainingSetClass &ts)
{

  // New models go here...add to the list and loop over paramters //
  if(strcmp(ts.model(),"TaylorF2_PN3pt5") == 0)
    TF2_FullWaveform(model_eval,params,xQuad,1.0,3.5); //amp=1.0,PN=3.5
  #ifdef MODEL_LAL
  else if (strcmp(ts.model(),"PhenomP_plus") == 0)
    PhenP_Waveform(model_eval, xQuad, params, ts.model());
  else if (strcmp(ts.model(),"PhenomP_cross") == 0)
    PhenP_Waveform(model_eval, xQuad, params, ts.model());
  else if(strcmp(ts.model(),"PhenomP_hphp") == 0)
    PhenP_Waveform(model_eval, xQuad, params, ts.model());
  else if (strcmp(ts.model(),"PhenomP_hchc") == 0)
    PhenP_Waveform(model_eval, xQuad, params, ts.model());
  else if (strcmp(ts.model(),"PhenomP_hphc") == 0)
    PhenP_Waveform(model_eval, xQuad, params, ts.model());
  else if (strcmp(ts.model(),"PhenomP_all_parts") == 0)
    PhenP_Waveform_All_Parts(model_eval, xQuad, params);
  else if (strcmp(ts.model(),"SEOBNRv2_ROM_SingleSpin") == 0)
    SEOBNRv2_ROM_SingleSpin_Waveform(model_eval, xQuad, params);
  #endif
  else {
    std::cerr << "Approximant not supported!" << std::endl;
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

  fprintf(stdout,"Populating training set on proc %i...\n",rank);
  fprintf(stdout,"Using the model %s\n",ts.model());

  gsl_vector_complex *model_eval;
  double *params;
  int proc_ts_size;

  params     = new double[ts.param_dim()]; // TF2 model param (mass 1, mass 2)
  model_eval = gsl_vector_complex_alloc(xQuad->size);

  ts.LocalTrainingSetSize(proc_ts_size,rank);

  for(int i = 0; i < proc_ts_size; i++){
    ts.GetParameterValue(params,rank,i); //params [global_i][j] * (scale[j])
    EvaluateModel(model_eval,xQuad,params,ts);
    gsl_matrix_complex_set_row(TS_gsl,i,model_eval);
  }

  // -- Normalize training space here -- //
  //fprintf(stdout,"Normalizing training set...\n");
  mygsl::NormalizeMatrixRows(TS_gsl,wQuad);

  delete[] params;
  gsl_vector_complex_free(model_eval);
}

};

#endif /* my_models.h */
