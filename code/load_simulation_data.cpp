#include <iostream>

#include "load_simulation_data.hpp"

#include "parameters.hpp"
#include "utils.h"
#include "gsl_helper_functions.hpp"
#include "quadratures.hpp"

LoadData::~LoadData() {

  if(load_eim_) {
    gsl_matrix_complex_free(eim_invV_);
    delete[] eim_indices_;
    delete[] eim_nodes_;
  }

  gsl_matrix_complex_free(RB_space_);
  gsl_vector_complex_free(wQuad_);
  gsl_vector_free(xQuad_);
  delete params_from_file_;
  params_from_file_ = NULL;

}

LoadData::LoadData(char** argv, bool load_eim) : 
  load_eim_(load_eim),
  base_path_(argv[2])
{

  params_from_file_ = new Parameters(argv,true);
  std::cout << *params_from_file_;

  // -- Deduce quadrature and basis sizes -- //
  rb_size_   = fcount_pts(argv[2],"/GreedyPoints.txt");
  if(params_from_file_->quad_type()==0 || params_from_file_->quad_type()==1){
    quad_size_ = params_from_file_->quad_points();
  } else {
    quad_size_ = fcount_pts(argv[2],"/quad_rule.txt");
  }

  if(load_eim_) {

    eim_invV_    = gsl_matrix_complex_alloc(rb_size_,rb_size_);
    eim_indices_ = new int[rb_size_];
    eim_nodes_   = new double[rb_size_];

    // -- Load EIM indices -- //
    std::string eim_ind_path(base_path_);
    eim_ind_path.append("EIM_indices.txt");
    FILE *fp_ind = fopen(eim_ind_path.c_str(), "r");
    if (fp_ind==NULL) {
      fprintf(stderr,"Could not open eim indices file.\n");
      exit(1);
    }
    int counter = 0;
    while(fscanf(fp_ind, "%i", &eim_indices_[counter]) != EOF)
      counter += 1;
    fclose(fp_ind);

    // -- Load EIM nodes -- //
    std::string eim_node_path(base_path_);
    eim_node_path.append("EIM_nodes.txt");
    FILE *fp_nodes = fopen(eim_node_path.c_str(), "r");
    if (fp_nodes==NULL) {
      fprintf(stderr,"Could not open eim idices file.\n");
      exit(1);
    }
    counter = 0;
    while(fscanf(fp_nodes, "%lf", &eim_nodes_[counter]) != EOF)
      counter += 1;

    fclose(fp_nodes);

    // -- Load EIM inverse Vandermonde -- //
    std::string eim_invV_path(base_path_);

    if (strcmp(argv[3],"gsl") == 0) {
      eim_invV_path.append("invV.bin");

      FILE *pf_invV = fopen(eim_invV_path.c_str(),"rb");
      if (pf_invV==NULL) {
        std::cerr << "could not open invV file\n";
        exit(1);
      }
      gsl_matrix_complex_fread(pf_invV,eim_invV_);
      fclose(pf_invV);
    }
    else if (strcmp(argv[3],"npy") == 0) {
      eim_invV_path.append("invV.npy");
      mygsl::gsl_matrix_complex_npy_load(eim_invV_path.c_str(),eim_invV_);
    }

  }

  RB_space_ = gsl_matrix_complex_alloc(rb_size_,quad_size_);

  // -- Read in basis from a binary file -- //
  std::string rb_filename(base_path_);
  std::cout << "reading basis from file...\n";

  if (strcmp(argv[3],"gsl") == 0) {

    rb_filename.append("/Basis.bin");
    FILE *pBASIS = fopen(rb_filename.c_str(),"rb");
    if (pBASIS==NULL) {
      std::cerr << "could not open basis file\n";
      exit(1);
    }
    gsl_matrix_complex_fread(pBASIS,RB_space_);
    fclose(pBASIS);

  }
  else if (strcmp(argv[3],"npy") == 0) {
    rb_filename.append("/Basis.npy");
    mygsl::gsl_matrix_complex_npy_load(rb_filename.c_str(),RB_space_);
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
