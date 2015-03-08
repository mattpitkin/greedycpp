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

    // overload constructors
    TrainingSetClass(Parameters *, int);
    TrainingSetClass(Parameters *, std::string); // random sampling. samples in file 

    // explicit destructor - constructor allocates memory on the heap
    ~TrainingSetClass();

    // BuildTS will decide how to populte the collection of points depending on input
    void BuildTS(const int *, const double *, const double *);
    void BuildTS(const char *);

    // Find the rank containing a training element via the global index //
    int FindRowIndxRank(const int) const;

    // Return local (on this proc) training set indexes //
    void GetLocalTrainingSet(int &, int &, int const);
    void LocalTrainingSetSize(int &, const int) const;

    // Return training set value (an array of size param_dim) //
    void GetParameterValue(double *,const int,const int) const;

    // accessor functions
    inline int ts_size() const { return ts_size_; }
    inline int param_dim() const { return param_dim_; }
    inline bool distributed() const { return distributed_; }
    inline const char * model() const { return model_; }
    inline const int * mystart() const {return mystart_; }
    inline const int * myend() const { return myend_; }
    inline const int * matrix_sub_size() const { return matrix_sub_size_; }
    inline const double * const* params() const { return params_; } // pointer to pointer also needs to be const
    inline const double * param_scale() const { return param_scale_; }

    // this routine writes ts to file //
    void WriteTrainingSet();

    // output ith (param_dim-dimensional) training set point //
    void fprintf_ith(FILE*,int) const;

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


    // allocates memory for params_ only called by construtors //
    void AllocTS();

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

