#ifndef LOADDATA_hpp
#define LOADDATA_hpp

#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_vector_complex_double.h>
class Parameters;

class LoadData {
  public:

    LoadData();
    LoadData(char **);
    ~LoadData();

    // accessor functions
    inline const
    Parameters& params_from_file() const { return *params_from_file_; }
    inline const int rb_size() const { return rb_size_; }
    inline const int quad_size() const { return quad_size_; }
    inline const gsl_matrix_complex& RB_space() const { return *RB_space_; }
    inline gsl_matrix_complex& RB_space() { return *RB_space_; } 
    inline const gsl_vector_complex& wQuad() const { return *wQuad_; }
    inline const gsl_vector& xQuad() const { return *xQuad_; }

  private:
    Parameters *params_from_file_;
    int rb_size_;
    int quad_size_;
    gsl_matrix_complex *RB_space_;
    gsl_vector_complex *wQuad_;
    gsl_vector *xQuad_;
};

#endif /* LoadData.hpp */
