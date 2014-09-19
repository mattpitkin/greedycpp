#ifndef training_set_hpp
#define training_set_hpp

#include "parameters.hpp"

/*
  Class contains all information about the training set and, if 
  MPI enabled, how its distributed over procs. 
*/


class TrainingSetClass {

    public:
        // default constructor -- not explicitly defined in cpp
        TrainingSetClass();

        // overload constructor
        TrainingSetClass(Parameters *, int);

        // explicit destructor - constructor allocates memory on the heap
        ~TrainingSetClass();

        // BuildTS will decide how to populte the collection of points depending on input
        void BuildTS(const int *, const double *, const double *);
        void BuildTS(const char *);

        // Description here TODO //
        int FindRowIndxRank(const int);

        // Return local (on this proc) training set indexes //
        void GetLocalTrainingSet(int &, int &, int const);
        void LocalTrainingSetSize(int &,const int);

        // Return training set value (an array of size param_dim) //
        void GetParameterValue(double *,const int,const int);

        // accessor functions
        inline int ts_size(){ return ts_size_; }
        inline int param_dim(){ return param_dim_; }
        inline bool distributed(){ return distributed_; }
        inline const char * model(){ return model_; }
        inline const int * mystart(){return mystart_; }
        inline const int * myend(){ return myend_; }
        inline const int * matrix_sub_size(){ return matrix_sub_size_; }
        inline const double * const* params(){ return params_; } // pointer to pointer also needs to be const
        inline const double * param_scale(){ return param_scale_; }

        // this routine writes ts to file //
        void WriteTrainingSet();

    private:

        // member variables
        double **params_;            // matrix where the columns are the parameters (param_dim columns) and the rows their indexing (ts_size rows)
        const double *param_scale_;  // see parameters.hpp (note: param_scale_ is actually a pointer to param_scale_ defined in Parameters class)
        int param_dim_;              // see parameters.hpp
        int ts_size_;                // see parameters.hpp
        char model_[100];            // see parameters.hpp
        bool distributed_;           // set to true if training space distributed over procs/nodes
        int *mystart_, *myend_;      // maps global row index of A onto local worker index 
        int *matrix_sub_size_;       // = ts.myend[rank]-ts.mystart[rank] = TS on each proc

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

