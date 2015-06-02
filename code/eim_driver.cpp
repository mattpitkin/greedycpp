// PURPOSE: Main executable for finding EIM points
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
  std::cout << "basis folder is: " << argv[2] << std::endl;
  std::cout << "basis format is: " << argv[3] << std::endl;

  std::string basis_path         = std::string(argv[2]);
  std::string random_sample_file = std::string(argv[3]);

  if(basis_path.at(basis_path.length()-1) != '/') {
    std::cerr << "directory with basis file must end with '/' " << std::endl;
    exit(1);
  }

  //--- Load greedy data from previous simulation -- //
  LoadData *data = new LoadData(argv);

  // -- build EIM information -- //
  EIM *eim = new EIM(data->rb_size(),data->quad_size(),false);
  eim->build_eim_data(&data->RB_space());

  delete data;
  data = NULL;

  LoadData *data_reload = new LoadData(argv);

  eim->rebuild_vandermonde(&data_reload->RB_space());
  eim->compute_vandermonde_inverse();

  const gsl_matrix_complex inv_vandermonde = eim->invV();
  const int *eim_indicies = eim->p();

  // add to output folder: relevant matrix, points (physical ones too?)
  //const gsl_vector xQuad = data->xQuad();

  // --- save data to file --- //
  // ...eim indicies
  std::string eim_indices_path(basis_path);
  eim_indices_path.append("EIM_indices.txt");
  FILE *eim_indices_file = fopen(eim_indices_path.c_str(),"w");
  for(int i=0;i<eim->eim_size();++i)
    fprintf(eim_indices_file,"%i\n",eim_indicies[i]);
  fclose(eim_indices_file);

  // ...inverse Vandermonde matrix
  //const char* datatype = data_reload->params_from_file().output_data_format().c_str();
  //std::cout << datatype << std::endl;
  std::string invV_path(basis_path);

  // txt
  //invV_path.append("invV");
  //mygsl::gsl_matrix_complex_fprintf(invV_path.c_str(),&eim->invV());

  // gsl bin
  //invV_path.append("invV.bin");
  //FILE *invV_data = fopen(invV_path.c_str(),"wb");
  //gsl_matrix_complex_fwrite(invV_data,&eim->invV());
  //fclose(invV_data);

  // numpy bin
  invV_path.append("invV.npy");
  mygsl::gsl_matrix_complex_npy_save(invV_path.c_str(),&eim->invV());

  std::string V_path(basis_path);
  V_path.append("V.npy");
  mygsl::gsl_matrix_complex_npy_save(V_path.c_str(),&eim->V());

  delete data_reload;
  data_reload = NULL;
  delete eim;
  eim = NULL;
}



