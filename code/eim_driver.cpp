// PURPOSE: Main executable for finding EIM points and the
//          associated interpolation matrix
// DATE: 5/31/2015


#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "load_simulation_data.hpp"
#include "parameters.hpp"
#include "gsl_helper_functions.hpp"
#include "eim.hpp"

int main (int argc, char **argv) {


  //----- Checking the number of Variables passed to the Executable -----//
  if (argc != 4) {
    std::cerr << 
         "Argument: 1. location of a cfg and parameter file" << std::endl;
    std::cerr << 
         "Argument: 2. directory (ending with /) with basis+quadrature info\n";
    std::cerr << "Argument: 3. basis format (npy or gsl)" << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "parameter file is: " << argv[1] << std::endl;
  std::cout << "basis folder is: "   << argv[2] << std::endl;
  std::cout << "basis format is: "   << argv[3] << std::endl;

  std::string basis_path         = std::string(argv[2]);
  std::string random_sample_file = std::string(argv[3]);

  if(basis_path.at(basis_path.length()-1) != '/') {
    std::cerr << "directory with basis file must end with '/' " << std::endl;
    exit(1);
  }

  //--- Load greedy data from previous simulation -- //
  LoadData *data = new LoadData(argv,false);

  // -- build EIM information -- //
  EIM *eim = new EIM(data->rb_size(),data->quad_size(),false);
  eim->build_eim_data(&data->RB_space());

  // -- reload the data -- eim has modifed RB_space inplace -- //
  delete data;
  data = new LoadData(argv,false);

  // -- compute the interpolation matrix -- //
  eim->rebuild_vandermonde(&data->RB_space());
  eim->compute_vandermonde_inverse();

  // --- save data to file --- //
  // ...eim indices
  std::string eim_indices_path(basis_path);
  eim_indices_path.append("EIM_indices.txt");
  FILE *eim_indices_file = fopen(eim_indices_path.c_str(),"w");
  for(int i=0;i<eim->eim_size();++i)
    fprintf(eim_indices_file,"%i\n",eim->p()[i]);
  fclose(eim_indices_file);

  // ...eim nodes
  const gsl_vector xQuad = data->xQuad();
  std::string eim_nodes_path(basis_path);
  eim_nodes_path.append("EIM_nodes.txt");
  FILE *eim_nodes_file = fopen(eim_nodes_path.c_str(),"w");
  for(int i=0;i<eim->eim_size();++i)
    fprintf(eim_nodes_file,"%f\n", gsl_vector_get(&xQuad,eim->p()[i]) );
  fclose(eim_indices_file);

  // ...inverse Vandermonde matrix
  std::string invV_path_txt(basis_path);
  std::string invV_path_gsl(basis_path);
  std::string invV_path_npy(basis_path);

  // save to txt file
  invV_path_txt.append("invV");
  mygsl::gsl_matrix_complex_fprintf(invV_path_txt.c_str(),&eim->invV());

  // save to gsl binary file 
  invV_path_gsl.append("invV.bin");
  FILE *invV_data = fopen(invV_path_gsl.c_str(),"wb");
  gsl_matrix_complex_fwrite(invV_data,&eim->invV());
  fclose(invV_data);

  // save to numpy binary file
  invV_path_npy.append("invV.npy");
  mygsl::gsl_matrix_complex_npy_save(invV_path_npy.c_str(),&eim->invV());

  delete data;
  data = NULL;
  delete eim;
  eim = NULL;
}



