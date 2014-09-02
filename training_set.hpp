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

// generates uniform spacing //
void uniform(const int &, const double &, const double &, double *);

// wrapper for 2D and ND tensor product builds //
void BuildTS_TensorProduct(const int *, const double *, const double *, TrainSet &);

// these routines will populate params in TrainSet -- 2-D //
void BuildTS_TensorProduct2D(const int *, const double *, const double *, TrainSet &);

// these routines will populate params in TrainSet -- N-D //
void BuildTS_RecursiveSetBuild(const double *params_low, const double *params_step_size, const int *params_num, const int level, double *param_vector, int &counter);
void BuildTS_TensorProductND(const int *params_num, const double *params_low, const double *params_high, TrainSet &ts);

// this routine writes ts to file //
void WriteTrainingSet(const TrainSet);


// --------------- work towards making this a class ------------------ //

class TrainingSetClass {

    public:
        // default constructor -- not explicitly defined in cpp
        TrainingSetClass();

        // overload constructor
        TrainingSetClass(int, double *, int, const char *);

        // destructor -- is default one ok if malloc called? what about new?
        //~TrainingSetClass();

        // BuildTS will decide how to populte the collection of points depending on input
        void BuildTS(int *, double *, double *);
        void BuildTS(const char *);

        // this routine distributed the ts over procs //
        void SplitTrainingSet(const int);

        // Descriptiong here TODO //
        int FindRowIndxRank(const int);

        //TrainSet GetTrainingSet() const;

        // accessor functions
        int ts_size();
        int param_dim();
        const char * model();
        bool distributed();
        const int* mystart();
        const int* myend();
        const int* matrix_sub_size();
        double** params(); // TODO: const
        const double* param_scale();

    private:
        // member variables
        //TrainSet ts_;

        double **params_;      // matrix where the columns are the parameters (param_dim columns) and the rows their indexing (ts_size rows)
        double *param_scale_;  // param_scale[i] scales ith paramter so that input to model is param_scale[i] * param_i (all scale's default is 1)
        int param_dim_;        // number of paramteric dimensions
        int ts_size_;          // number of training set elements
        char model_[100];      // name of model (ex: TaylorF2_PN3pt5)

        // settings for MPI runs //
        bool distributed_;     // set to true if TS distributed over procs/nodes
        int *mystart_, *myend_; // maps global row index of A onto local worker index 
        int *matrix_sub_size_; // = ts.myend[rank]-ts.mystart[rank] = TS on each proc


        // wrapper for 2D and ND tensor product builds //
        //void BuildTS_TensorProduct2D(const int *, const double *, const double *);

        // these routines will populate params in TrainSet -- 2-D //
        //void BuildTS_TensorProduct2D(const int *, const double *, const double *);
        void BuildTS_FromFile2D(const char *);
        void BuildTS_FromFileND(const char *);

        // these routines will populate params in TrainSet -- N-D //
        //void BuildTS_RecursiveSetBuild(const double *, const double *, const int *, const int, double *, int &);
        //void BuildTS_TensorProductND(const int *, const double *, const double *);


};

#endif /* training_set_hpp */

