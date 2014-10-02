// AUTHOR :  Scott Field 
//           sfield@astro.cornell.edu
//
// DATE: Oct 1, 2014
//
// PURPOSE: test basis accuracy from random model samples

//#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <cmath>

#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_complex.h>
#include <gsl/gsl_block_complex_float.h>

#include "parameters.hpp"
#include "gsl_helper_functions.h"

int main (int argc, char **argv) {

  //----- Checking the number of Variables passed to the Executable -----//
  if (argc != 2) {
    std::cerr << "Argument: 1. location of a cfg configuration/parameter file (ends in .cfg)" << std::endl;
    std::cerr << "Argument: 2. directory containing basis and quadrature information" << std::endl;
    //exit(0);
  }
  std::cout << "parameter file is: " << argv[1] << std::endl;
  std::cout << "basis folder is: " << argv[2] << std::endl;

  //--- Read input file argv[1]. If there is an error, report and exit.
  //--- Parameters class contains relevant information about parameters 
  Parameters *params_from_file = new Parameters(argv);
  std::cout << *params_from_file;

  gsl_matrix_complex *RB_space;
  RB_space = gsl_matrix_complex_alloc(params_from_file->max_RB(),params_from_file->quad_points());

  // create basis verifcation folder as subfolder of output

  // read in basis and output to file
  char rb_filename[100];
  FILE *pBASIS;
  strcpy(rb_filename,argv[2]);
  strcat(rb_filename,"/Basis.bin");

  std::cout << "reading basis from file " << rb_filename << " ..." << std::endl;

  pBASIS = fopen(rb_filename,"rb");
  gsl_matrix_complex_fread(pBASIS,RB_space);
  fclose(pBASIS);

  //char testwrite[100];
  //strcpy(testwrite,"RB_test");
  //mygsl::gsl_matrix_complex_fprintf(testwrite,RB_space);



  gsl_matrix_complex_free(RB_space);

  delete params_from_file;
  params_from_file = NULL;
}

