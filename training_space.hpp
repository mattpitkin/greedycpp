#ifndef training_space_hpp
#define training_space_hpp

#include "training_set.hpp"
#include "parameters.hpp"

/*
  Class contains all information about the training space.
*/


class TrainingSpaceClass : public TrainingSetClass {
public:

    TrainingSpaceClass(Parameters *,int);

    //~TrainingSpaceClass();

};

#endif /* training_set_hpp */

