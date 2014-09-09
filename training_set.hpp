#ifndef training_set_hpp
#define training_set_hpp

/*
  Class contains all information about the training set and, if 
  MPI enabled, how its distributed over procs. 
*/


class TrainingSetClass {

    public:
        // default constructor -- not explicitly defined in cpp
        TrainingSetClass();

        // overload constructor
        TrainingSetClass(int, double *, int, const char *,int);

        // destructor -- TODO: is default one ok if malloc called? what about new?
        //~TrainingSetClass();

        // BuildTS will decide how to populte the collection of points depending on input
        void BuildTS(int *, double *, double *);
        void BuildTS(const char *);

        // Description here TODO //
        int FindRowIndxRank(const int);

        // Return local (on this proc) training set indexes //
        void GetLocalTrainingSet(int &, int &, int const);
        void LocalTrainingSetSize(int &,const int);

        // Return training set value (an array of size param_dim) //
        void GetParameterValue(double *,const int,const int);

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

        // this routine writes ts to file //
        void WriteTrainingSet();

    private:

        // member variables
        double **params_;       // matrix where the columns are the parameters (param_dim columns) and the rows their indexing (ts_size rows)
        double *param_scale_;   // param_scale[i] scales ith paramter so that input to model is param_scale[i] * param_i (all scale's default is 1)
        int param_dim_;         // number of paramteric dimensions
        int ts_size_;           // number of training set elements
        char model_[100];       // name of model (ex: TaylorF2_PN3pt5)
        bool distributed_;      // set to true if TS distributed over procs/nodes
        int *mystart_, *myend_; // maps global row index of A onto local worker index 
        int *matrix_sub_size_;  // = ts.myend[rank]-ts.mystart[rank] = TS on each proc

        // this routine distributed the ts over procs //
        void SplitTrainingSet(const int);

        // generates uniform spacing //
        void uniform(const int &, const double &, const double &, double *);

        // these routines will populate params in TrainSet -- 2-D //
        void BuildTS_TensorProduct2D(const int *, const double *, const double *);

        // wrapper for 2D and ND tensor product builds //
        void BuildTS_TensorProduct(const int *, const double *, const double *);

        // these routines will populate params in TrainSet -- N-D //
        void BuildTS_RecursiveSetBuild(const double *, const double *, const int *, const int, double *, int &);
        void BuildTS_TensorProductND(const int *, const double *, const double *);

        // these routines will populate params in TrainSet -- 2-D //
        void BuildTS_FromFile2D(const char *);
        void BuildTS_FromFileND(const char *);
};

#endif /* training_set_hpp */

