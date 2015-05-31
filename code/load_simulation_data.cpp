#include <iostream>

#include "load_simulation_data.hpp"

#include "parameters.hpp"
#include "utils.h"
#include "gsl_helper_functions.hpp"
#include "quadratures.hpp"

LoadData::~LoadData() {

  gsl_matrix_complex_free(RB_space_);
  gsl_vector_complex_free(wQuad_);
  gsl_vector_free(xQuad_);
  delete params_from_file_;
  params_from_file_ = NULL;

}

LoadData::LoadData(char** argv) {

  params_from_file_ = new Parameters(argv,true);
  std::cout << *params_from_file_;

  // -- Deduce quadrature and basis sizes -- //
  //int rb_size   = fcount_pts(argv[2],"/ApproxErrors.txt");
  rb_size_   = fcount_pts(argv[2],"/GreedyPoints.txt");
  if(params_from_file_->quad_type()==0 || params_from_file_->quad_type()==1){
    quad_size_ = params_from_file_->quad_points();
  } else {
    quad_size_ = fcount_pts(argv[2],"/quad_rule.txt");
  }

  RB_space_ = gsl_matrix_complex_alloc(rb_size_,quad_size_);

  // -- Read in basis from a binary file -- //
  char rb_filename[200];
  std::cout << "reading basis from file...\n";
  strcpy(rb_filename,argv[2]);

  if (strcmp(argv[4],"gsl") == 0) {

    strcat(rb_filename,"/Basis.bin");
    FILE *pBASIS;
    pBASIS = fopen(rb_filename,"rb");
    if (pBASIS==NULL) {
      std::cerr << "could not open basis file\n";
      exit(1);
    }
    gsl_matrix_complex_fread(pBASIS,RB_space_);
    fclose(pBASIS);

  }
  else if (strcmp(argv[4],"npy") == 0) {
    strcat(rb_filename,"/Basis.npy");
    mygsl::gsl_matrix_complex_npy_load(rb_filename,RB_space_);
  }
  else {
    std::cerr << "file type not supported\n";
    exit(1);
  }

  // -- reconstruct quadrature rule used in greedycpp -- //
  SetupQuadratureRule(&wQuad_,&xQuad_,params_from_file_);
  assert(quad_size_==xQuad_->size);
  assert(quad_size_==params_from_file_->quad_points());

}
