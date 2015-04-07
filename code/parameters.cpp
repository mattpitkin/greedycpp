#include <cstdlib>
#include <iostream>
#include <libconfig.h++>
#include <fstream>
#include <string.h>
#include <cmath>
//#include <time.h>
//#include <vector>
#include <sstream>

#include "parameters.hpp"
#include "utils.h"


Parameters::~Parameters(){

  // TODO: this gets deleted inside trainingset class. should copy pointer to ts class instead
  //delete[] param_scale_;

  if(!load_from_file_) {
    delete[] params_num_;
    delete[] params_low_;
    delete[] params_high_;
  }

}

Parameters::Parameters(char ** argv, bool high_verbosity){

  //--- Read input (config) file. If there is an error, report and exit ---//
  // TODO: with MPI, multiple procs reading same paramter file... seems bad //
  libconfig::Config cfg;
  try{
    cfg.readFile(argv[1]);
  }
  catch(const libconfig::FileIOException &fioex){
    std::cerr << "Error while reading cfg file. Check file's name and path."
              << std::endl;
    exit(1);
  }
  catch(const libconfig::ParseException &pex){
    std::cerr << "Parse error. Check cfg file for correct syntax." 
              << std::endl;
    exit(1);
  }

  bool cfg_status;
  config_file_name_ == std::string(argv[1]); 
  const libconfig::Setting& root = cfg.getRoot();

  if(high_verbosity) {
    fprintf(stdout,"Loading seed, tol, max_RB, weighted_inner, param_dim, "
                    "export_ts, load_from_file, quad_type...\n");
  }
  load_from_file_ = cfg.lookup("load_from_file");
  quad_type_      = cfg.lookup("quad_type");
  seed_           = cfg.lookup("seed");
  tol_            = cfg.lookup("tol");
  max_RB_         = cfg.lookup("max_RB");
  weighted_inner_ = cfg.lookup("weighted_inner");
  param_dim_      = cfg.lookup("param_dim");
  export_tspace_  = cfg.lookup("export_tspace");

  if(high_verbosity) {
    fprintf(stdout,"Loading param_scale...\n");
  }
  param_scale_ = new double[param_dim_];
  char scale_str[20];
  for(int i = 0; i < param_dim_; i++){
    snprintf(scale_str, 20, "p%d_scale", i+1);
    param_scale_[i] =  cfg.lookup(scale_str);
  }

  if(high_verbosity) {
    fprintf(stdout,"Loading model_name, output_dir, output_data_format...\n");
  }
  if(cfg.lookupValue("model_name",model_name_)
   && cfg.lookupValue("output_dir",output_dir_)
   && cfg.lookupValue("output_data_format",output_data_format_)){
    if(high_verbosity) {
      fprintf(stdout,"Successfully loaded model name, output directory "
                     "location and output data format type\n");
      }
  }
  else{
    fprintf(stderr,"Failed to load either model name, output directory "
                    "location or output data format type\n");
    exit(1);
  }


  // -- build training set (builds list of paramter values in ts.params) -- //
  if (load_from_file_) {
    if(high_verbosity) {
      fprintf(stdout,"Loading ts values from file...\n");
    }
    cfg_status = cfg.lookupValue("ts_file", ts_file_name_);
    if (!cfg_status){
      fprintf(stderr, "ts_file not found in config file\n");
      exit(1);
    }
    ts_size_ = fcount_pts(ts_file_name_.c_str()); // assumes each row is point
    if(high_verbosity) { 
      std::cout << "training set file found to " 
                << ts_size_ << " parameter samples" << std::endl;
    }
  }
  else {
    if(high_verbosity) {
      fprintf(stdout,"Loading param ranges from file...\n");
    }
    params_num_       = new int[param_dim_];
    params_low_       = new double[param_dim_];
    params_high_      = new double[param_dim_];

    // read in params_num and determine ts_size //
    ts_size_ = 1;
    libconfig::Setting& params_num_s = root["params_num"];
    int count_array = params_num_s.getLength();
    if(count_array != param_dim_){
      fprintf(stderr,"elements in params_num defined in configuration file does equal param_dim\n");
      exit(1);
    }
    for(int i = 0; i < param_dim_; i++){
      params_num_[i] = params_num_s[i];
      ts_size_       = params_num_[i]*ts_size_; // assumes tensor product TS
    }

    // read in params_low //
    libconfig::Setting& params_low_s = root["params_low"];
    count_array = params_low_s.getLength();
    if(count_array != param_dim_){
      fprintf(stderr,"elements in params_low defined in configuration file does equal param_dim\n");
      exit(1);
    }
    for(int i = 0; i < param_dim_; i++){
      params_low_[i] = params_low_s[i];
    }

    // read in params_high //
    libconfig::Setting& params_high_s = root["params_high"];
    count_array = params_high_s.getLength();
    if(count_array != param_dim_){
      fprintf(stderr,"elements in params_high defined in configuration file does equal param_dim\n");
      exit(1);
    }
    for(int i = 0; i < param_dim_; i++){
      params_high_[i] = params_high_s[i];
    }
  }

  if(high_verbosity) {
    fprintf(stdout,"Reading quadrature information from cfg file...\n");
  }
  if (quad_type_ == 2) // user specified quadrature points, load these from file
  {
    cfg_status = cfg.lookupValue("quad_nodes_file", quad_nodes_file_);
    if (!cfg_status){
      fprintf(stderr, "quadrature node file name must be specified\n");
      exit(1);
    }

    quad_points_ = fcount_pts(quad_nodes_file_.c_str());

    cfg_status = cfg.lookupValue("num_weight_file", num_weight_file_);
    if (!cfg_status){
      fprintf(stderr, "quadrature weight file name must be specified\n");
       exit(1);
    }
  }
  else
  {
    x_min_       = cfg.lookup("x_min");
    x_max_       = cfg.lookup("x_max");
    quad_points_ = cfg.lookup("quad_points");
  }

  if (weighted_inner_){
    cfg_status = cfg.lookupValue("quad_weight_file", quad_weight_file_);
    if (!cfg_status){
      fprintf(stderr, "weight_file not found in config file\n");
      exit(1);
    }
  }

}

std::ostream& operator<<(std::ostream& os, const Parameters& p)
{
  os << "Dimensionality of parameter space " 
     << p.param_dim() << " sampled by " << p.ts_size() 
     << " training set points.\n"
     << "Using waveform model " << p.model_name() << std::endl;
  return os;
}


// Note: inlined accessor methods located in hpp file //
