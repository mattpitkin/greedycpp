#ifndef TrainingSet_hpp
#define TrainingSet_hpp

/*
  TrainSet should contain all information about the training set and, if 
  MPI enabled, how its distributed over procs. 
*/
typedef struct
tagTrainSet
{
  double *m1, *m2;      // each element of the training set is (m1[i],m2[i])
  double p1_scale;      // scale each m1[i] so that input to model is p1_scale * m1[i] (default p1_scale = 1)
  double p2_scale;      // same as p1_scale, but for m2[i]
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

// these routines will populate m1 and m2 //
void BuildTS_tensor_product(const int &, const double &, const double &, TrainSet &);
void BuildTS_from_file(const char *, TrainSet &);

// this routine distributed the ts over procs //
void SplitTrainingSet(const int, TrainSet &);

// this routine writes ts to file //
void WriteTrainingSet(const TrainSet);


#endif /* TrainingSet_hpp */

