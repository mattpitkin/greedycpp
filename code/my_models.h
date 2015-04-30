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
#include "gsl_helper_functions.h"

// model specific header files //
#include "../models/taylorf2/spa_waveforms.h"
#include "../models/lal/phenomp.h"


namespace mymodel {

void EvaluateModel(gsl_vector_complex *model_eval,
                     const gsl_vector *xQuad,
                     const double *params,
                     const TrainingSetClass &ts)
{


  // *** BEGIN MODEL SPECIFIC SECTION *** //
  // New models go here...add to the list and loop over paramters //
  if(strcmp(ts.model(),"TaylorF2_PN3pt5") == 0){
    //fprintf(stdout,"Using the TaylorF2 spa approximant to PN=3.5\n");
    TF2_FullWaveform(model_eval,params,xQuad,1.0,3.5); //amp=1.0,PN=3.5
  }
  else if (strcmp(ts.model(),"PhenomP_plus") == 0){
    //fprintf(stdout,"Using the PhenP approx\n");
    PhenP_Waveform(model_eval, xQuad, params, ts.model());
  }
  else if (strcmp(ts.model(),"PhenomP_cross") == 0){
    //fprintf(stdout,"Using the PhenP approx\n");
    PhenP_Waveform(model_eval, xQuad, params, ts.model());
  }
  else if(strcmp(ts.model(),"PhenomP_hphp") == 0){
    //fprintf(stdout,"Using the PhenP  approx (hplus squared)\n");
    PhenP_Waveform(model_eval, xQuad, params, ts.model());
  }
  else if (strcmp(ts.model(),"PhenomP_hchc") == 0){
    //fprintf(stdout,"Using the PhenP  approx (hcross squared)\n");
    PhenP_Waveform(model_eval, xQuad, params, ts.model());
  }
  else if (strcmp(ts.model(),"PhenomP_hphc") == 0){
    //fprintf(stdout,"Using the PhenP  approx (hplus x hcross)\n");
    PhenP_Waveform(model_eval, xQuad, params, ts.model());
  }
  else{
    std::cerr << "Approximant not supported!" << std::endl;
    exit(1);
  }
}

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

  // *** BEGIN MODEL SPECIFIC SECTION *** //
  for(int i = 0; i < proc_ts_size; i++){
    ts.GetParameterValue(params,rank,i); //params [global_i][j] * (scale[j])
    EvaluateModel(model_eval,xQuad,params,ts);
    gsl_matrix_complex_set_row(TS_gsl,i,model_eval);
  }
  /*// New models go here...add to the list and loop over paramters //
  if(strcmp(ts.model(),"TaylorF2_PN3pt5") == 0){
    fprintf(stdout,"Using the TaylorF2 spa approximant to PN=3.5\n");

    for(int i = 0; i < proc_ts_size; i++){
      ts.GetParameterValue(params,rank,i); //params [global_i][j] * (scale[j])
      TF2_FullWaveform(model_eval,params,xQuad,1.0,3.5); //amp=1.0,PN=3.5
      gsl_matrix_complex_set_row(TS_gsl,i,model_eval);
    }

  }
  else if (strcmp(ts.model(),"PhenomP_plus") == 0){
    fprintf(stdout,"Using the PhenP approx\n");

    for(int i = 0; i < proc_ts_size; i++){
      ts.GetParameterValue(params,rank,i);
      PhenP_Waveform(model_eval, xQuad, params, ts.model());
      gsl_matrix_complex_set_row(TS_gsl,i,model_eval);
    }
  }
  else if (strcmp(ts.model(),"PhenomP_cross") == 0){
    fprintf(stdout,"Using the PhenP approx\n");

    for(int i = 0; i < proc_ts_size; i++){
      ts.GetParameterValue(params,rank,i);
      PhenP_Waveform(model_eval, xQuad, params, ts.model());
      gsl_matrix_complex_set_row(TS_gsl,i,model_eval);
    }
  }
  else if(strcmp(ts.model(),"PhenomP_hphp") == 0){
    fprintf(stdout,"Using the PhenP  approx (hplus squared)\n");

    for(int i = 0; i < proc_ts_size; i++){
      ts.GetParameterValue(params,rank,i);
      PhenP_Waveform(model_eval, xQuad, params, ts.model());
      gsl_matrix_complex_set_row(TS_gsl,i,model_eval);
    }
  }  
  else if (strcmp(ts.model(),"PhenomP_hchc") == 0){
    fprintf(stdout,"Using the PhenP  approx (hcross squared)\n");

    for(int i = 0; i < proc_ts_size; i++){
      ts.GetParameterValue(params,rank,i);
      PhenP_Waveform(model_eval, xQuad, params, ts.model());
      gsl_matrix_complex_set_row(TS_gsl,i,model_eval);
    }
  }
  else if (strcmp(ts.model(),"PhenomP_hphc") == 0){
    fprintf(stdout,"Using the PhenP  approx (hplus x hcross)\n");

    for(int i = 0; i < proc_ts_size; i++){
      ts.GetParameterValue(params,rank,i);
      PhenP_Waveform(model_eval, xQuad, params, ts.model());
      gsl_matrix_complex_set_row(TS_gsl,i,model_eval);
    }
  }
  else{
    std::cerr << "Approximant not supported!" << std::endl;
    exit(1);
  }*/
  // *** END MODEL SPECIFIC SECTION *** //


  // -- Normalize training space here -- //
  //fprintf(stdout,"Normalizing training set...\n");
  mygsl::NormalizeMatrixRows(TS_gsl,wQuad);

  delete[] params;
  gsl_vector_complex_free(model_eval);
}

};

#endif /* my_models.h */
