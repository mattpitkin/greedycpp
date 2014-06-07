#ifndef TrainingSet_H
#define TrainingSet_H

/*

  TrainSet data structure is used by greedy algorithm. TrainSet should
  contain all information about the training set and, if using MPI,
  how its distributed over procs. 

*/


/*-- Training Set data structure -- */
typedef struct
tagTrainSet
{
  double *m1, *m2;      // each element of the training set is (m1[i],m2[i])
  int param_dim;        // number of paramteric dimensions
  int ts_size;          // number of training set elements
  bool distributed;     // true if training set will be split over worker procs/nodes
  int *matrix_sub_size; // ts.myend[rank] - ts.mystart[rank] --> number of rows on each proc
  int *mystart, *myend; // maps global row index of A into local worker index (when distributed = true)
  char model[100];      // name of model (ex: TaylorF2_PN3pt5)
}
TrainSet;

// determine training set sizes and allocates memory //
void ts_alloc(const int ts_size, const int param_dim, const char *model_name, TrainSet &ts);

// generates uniform spacing //
void uniform(const int &, const double &, const double &, double *);

// these routines will populate m1 and m2 //
void BuildTS_tensor_product(const int &, const double &, const double &, TrainSet &);
void BuildTS_from_file(const char *, TrainSet &);

// this routine distributed the ts over procs //
void SplitTrainingSet(const int, TrainSet &);

// this routine writes ts to file //
void WriteTrainingSet(const TrainSet);


#endif

