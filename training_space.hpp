#ifndef training_space_hpp
#define training_space_hpp

#include "training_set.hpp"

/*
  Class contains all information about the training space.
*/


class TrainingSpaceClass : public TrainingSetClass {
public:

    TrainingSpaceClass(int, double *, int, const char *,int);

    //~TrainingSpaceClass();

};

#endif /* training_set_hpp */

