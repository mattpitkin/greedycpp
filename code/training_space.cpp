#include <iostream>
#include <fstream>
#include <string.h>
#include <cmath>
#include <time.h>
#include <vector>
#include <sstream>

#include "training_set.hpp"
#include "training_space.hpp"
#include "parameters.hpp"


/*
TrainingSpaceClass::~TrainingSpaceClass() {

    for(int j = 0; j < ts_size_; j++) {
      delete [] params_[j];
    }
    delete [] params_;

    delete [] param_scale_;
    delete [] mystart_;
    delete [] myend_;
    delete [] matrix_sub_size_;

}
*/

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

