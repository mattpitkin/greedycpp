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
#include "models/spa_waveforms.h"


namespace mymodel {

void FillTrainingSet(gsl_matrix_complex *TS_gsl,
                     const gsl_vector *xQuad,
                     const gsl_vector_complex *wQuad,
                     const TrainingSetClass &ts,
                     const int rank)
{

  fprintf(stdout,"Populating training set on proc %i...\n",rank);

  gsl_vector_complex *model_eval;
  double *params;
  int proc_ts_size;

  params     = new double[ts.param_dim()]; // TF2 model param (mass 1, mass 2)
  model_eval = gsl_vector_complex_alloc(xQuad->size);

  ts.LocalTrainingSetSize(proc_ts_size,rank);

  // *** BEGIN MODEL SPECIFIC SECTION *** //
  // New models go here...add to the list and loop over paramters //
  if(strcmp(ts.model(),"TaylorF2_PN3pt5") == 0){

    fprintf(stdout,"Using the TaylorF2 spa approximant to PN=3.5\n");

    for(int i = 0; i < proc_ts_size; i++){
      // fills params at [global_i][j] * (param_scale[j]) //
      ts.GetParameterValue(params,rank,i);
      TF2_FullWaveform(model_eval,params,xQuad,1.0,3.5); //amp=1.0,PN=3.5
      gsl_matrix_complex_set_row(TS_gsl,i,model_eval);
    }

  }
  else{
    std::cerr << "Approximant not supported!" << std::endl;
    exit(1);
  }
  // *** END MODEL SPECIFIC SECTION *** //


  // -- Normalize training space here -- //
  fprintf(stdout,"Normalizing training set...\n");
  mygsl::NormalizeMatrixRows(TS_gsl,wQuad);

  delete[] params;
  gsl_vector_complex_free(model_eval);
}

};

#endif /* my_models.h */
