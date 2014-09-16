#include <iostream>
#include <fstream>
#include <string.h>
#include <cmath>
#include <time.h>
#include <vector>
#include <sstream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_complex.h>
#include <gsl/gsl_block_complex_float.h>

#include "training_set.hpp"
#include "training_space.hpp"
#include "parameters.hpp"

TrainingSpaceClass::TrainingSpaceClass(Parameters *p, int procs_size)
    : TrainingSetClass(p,procs_size)
{
    // this could make more sense in training set. keep here for now
    // fills params_ with values //
    if(p->load_from_file()){
        BuildTS(p->ts_file_name().c_str());
    }
    else{
        BuildTS(p->params_num(),p->params_low(),p->params_high());
    }

 
}

