#ifndef LOADDATA_hpp
#define LOADDATA_hpp

#include <cassert>
#include <string>

#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_vector_complex_double.h>
class Parameters;

class LoadData {
  public:

    LoadData();
    LoadData(char ** argv,bool load_eim);
    ~LoadData();

    // accessor functions
    inline const
    Parameters& params_from_file() const { return *params_from_file_; }
    inline int rb_size() const { return rb_size_; }
    inline int quad_size() const { return quad_size_; }
    inline const gsl_matrix_complex& RB_space() const { return *RB_space_; }
    inline gsl_matrix_complex& RB_space() { return *RB_space_; } 
    inline const gsl_vector_complex& wQuad() const { return *wQuad_; }
    inline const gsl_vector& xQuad() const { return *xQuad_; }

    // accessor to eim data -- available if load_eim_ = true
    //inline const int eim_indices(int i) const { return eim_indices_[i]; }
    inline int* eim_indices() {return eim_indices_; }
    inline double eim_nodes(int i) { return eim_nodes_[i]; }
    inline gsl_matrix_complex* invV() { return eim_invV_; }


  private:
    Parameters *params_from_file_;
    int rb_size_;
    int quad_size_;
    gsl_matrix_complex *RB_space_;
    gsl_vector_complex *wQuad_;
    gsl_vector *xQuad_;
    int *eim_indices_;
    double *eim_nodes_;
    gsl_matrix_complex *eim_invV_;
    bool load_eim_;
    std::string base_path_;
};

#endif /* LoadData.hpp */
