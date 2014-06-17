#ifndef training_set_hpp
#define training_set_hpp

/*
  TrainSet should contain all information about the training set and, if 
  MPI enabled, how its distributed over procs. 
*/
typedef struct
tagTrainSet
{
  double **params;      // matrix where the columns are the parameters (param_dim columns) and the rows their indexing (ts_size rows)
  double *param_scale;  // param_scale[i] scales ith paramter so that input to model is param_scale[i] * param_i (all scale's default is 1)
  int param_dim;        // number of paramteric dimensions
  int ts_size;          // number of training set elements
  char model[100];      // name of model (ex: TaylorF2_PN3pt5)

  // settings for MPI runs //
  bool distributed;     // set to true if TS distributed over procs/nodes
  int *mystart, *myend; // maps global row index of A onto local worker index 
  int *matrix_sub_size; // = ts.myend[rank]-ts.mystart[rank] = TS on each proc 

}
TrainSet;

// determine training set size and allocate memory //
void ts_alloc(const int ts_size, const int param_dim, 
              const char *model_name, TrainSet &ts);

// generates uniform spacing //
void uniform(const int &, const double &, const double &, double *);

// these routines will populate params in TrainSet //
void BuildTS_TensorProduct2D(const int *, const double *, const double *, TrainSet &);
void BuildTS_FromFile(const char *, TrainSet &);

// this routine distributed the ts over procs //
void SplitTrainingSet(const int, TrainSet &);

// this routine writes ts to file //
void WriteTrainingSet(const TrainSet);


#endif /* training_set_hpp */

