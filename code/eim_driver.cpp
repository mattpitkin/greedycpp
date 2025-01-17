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
  std::cout << "Finished loading basis.\n";

  // -- compute the interpolation matrix -- //
  eim->rebuild_vandermonde(&data->RB_space());
  eim->compute_vandermonde_inverse();

  // --- save data to file --- //
  // ...eim indices (txt only)
  // ...eim nodes (txt only)
  clock_t start, end;
  start = clock();
  const gsl_vector xQuad = data->xQuad();
  std::cout << "Saving EIM indices and nodes...\n";
  std::string eim_indices_path(basis_path);
  std::string eim_nodes_path(basis_path);
  eim_indices_path.append("EIM_indices.txt");
  eim_nodes_path.append("EIM_nodes.txt");
  FILE *eim_indices_file = fopen(eim_indices_path.c_str(),"w");
  FILE *eim_nodes_file   = fopen(eim_nodes_path.c_str(),"w");
  for(int i=0;i<eim->eim_size();++i)
  {
    const int p_eim = eim->p()[i];
    fprintf(eim_indices_file,"%i\n",p_eim);
    fprintf(eim_nodes_file,"%f\n", gsl_vector_get(&xQuad,p_eim) );
  }
  fclose(eim_indices_file);
  fclose(eim_nodes_file);
  end = clock();
  double alg_time = ((double) (end - start)/CLOCKS_PER_SEC);
  fprintf(stdout,"Saving EIM indices and nodes took %f seconds\n",alg_time);


  // ...inverse Vandermonde matrix (txt, numpy or gsl)
  std::cout << "Saving EIM inverse interpolation matrix...\n";
  start = clock();  
  bool wrote = false;
  std::string datatype(data->params_from_file().output_data_format());

  if(datatype.compare("txt")==0 || datatype.compare("both")==0) {
    std::string invV_path_txt(basis_path);
    invV_path_txt.append("invV");
    mygsl::gsl_matrix_complex_fprintf(invV_path_txt.c_str(),&eim->invV());
    wrote = true;
  }
  if(datatype.compare("bin")==0 || datatype.compare("both")==0) {
    std::string invV_path_gsl(basis_path);
    invV_path_gsl.append("invV.bin");
    std::cout << "Saving file to " << invV_path_gsl << " ...\n";
    FILE *invV_data = fopen(invV_path_gsl.c_str(),"wb");
    if (invV_data==NULL) 
    {
      std::cout << "ERROR: could not open file\n";
    }
    gsl_matrix_complex_fwrite(invV_data,&eim->invV());
    fclose(invV_data);
    wrote = true;
  }
  if(datatype.compare("npy")==0){
    std::string invV_path_npy(basis_path);
    invV_path_npy.append("invV.npy");
    mygsl::gsl_matrix_complex_npy_save(invV_path_npy.c_str(),&eim->invV());
    wrote = true;
  }
  if(!wrote){
    fprintf(stderr,"file type not supported");
    exit(1);
  }
  end = clock();
  alg_time = ((double) (end - start)/CLOCKS_PER_SEC);
  fprintf(stdout,"Saving EIM matrix took %f seconds\n",alg_time);

  delete data;
  data = NULL;
  delete eim;
  eim = NULL;
}



