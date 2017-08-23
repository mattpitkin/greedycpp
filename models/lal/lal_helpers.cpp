// AUTHOR :  Scott Field 
//           sfield@astro.cornell.edu
//
// DATE: May 30, 2016
//
// PURPOSE: various lal helper routines

#include <cmath>
#include <iostream>
#include <cstdlib>

#include "lal_helpers.hpp"

// function to split a string on a delimiter (from http://stackoverflow.com/a/236803/1862861)
void split_string(const std::string &s, char delim, std::vector<std::string> &elems)
{
  std::stringstream ss;
  ss.str(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
}

std::vector<std::string> split_string(const std::string &s, char delim)
{
  std::vector<std::string> elems;
  split_string(s, delim, elems);
  return elems;
}

namespace lal_help {

// TODO: use enum instead of string tags

std::string get_waveform_part_tag(const std::string model_tag)
{
  //std::size_t pos  = model_tag.find("_");
  std::size_t pos  = model_tag.find_last_of("_");
  std::string part = model_tag.substr(pos+1);
  return part;
}

std::string get_waveform_part_tag(const double model_type)
{
  std::string model_part;
  if( std::abs( model_type - 0) < 1.e-10) {
    //fprintf(stdout,"Plus with 7th param %f\n",model_type);
    //PhenP_Waveform(wv,fnodes,params,"PhenomP_plus");
    model_part = "plus";
  }
  else if( std::abs( model_type - 1) < 1.e-10) {
    //fprintf(stdout,"cross with 7th param %f\n",model_type);
    //PhenP_Waveform(wv,fnodes,params,"PhenomP_cross");
    model_part = "cross";
  }
  else if( std::abs( model_type - 2) < 1.e-10) {
    //PhenP_Waveform(wv,fnodes,params,"PhenomP_hphp");
    model_part = "hphp";
  }
  else if( std::abs( model_type - 3) < 1.e-10) {
    //PhenP_Waveform(wv,fnodes,params,"PhenomP_hchc");
    model_part = "hchc";
  }
  else if( std::abs( model_type - 4) < 1.e-10) {
    //PhenP_Waveform(wv,fnodes,params,"PhenomP_hphc");
    model_part = "hphc";
  }
  else if( std::abs( model_type - 5) < 1.e-10) {
    model_part = "hpPLUShcSquared";
  }
  else {
    std::cerr << "PhenomP all parts -- unknown part" << std::endl;
    exit(1);
  }
  return model_part;
}


void lal_waveform_part(gsl_vector_complex *wv,
                       const std::string model_part,
                       const COMPLEX16FrequencySeries *hptilde,
                       const COMPLEX16FrequencySeries *hctilde,
                       const int n)
{

  if(model_part.compare("plus") == 0) {
    //fprintf(stdout,"Im plus\n");
    for (int i=0; i<n; i++) {
      gsl_vector_complex_set(wv, i, (hptilde->data->data)[i]);
    }
  }
  else if(model_part.compare("cross") == 0) {
    //fprintf(stdout,"Im cross\n");
    for (int i=0; i<n; i++) {
      gsl_vector_complex_set(wv, i, (hctilde->data->data)[i]);
    }
  }
  else if(model_part.compare("hphp") == 0) {
    //fprintf(stdout,"hp conj(hp)\n");
    for (int i=0; i<n; i++) {
      gsl_vector_complex_set(wv, i,
        gsl_complex_mul((hptilde->data->data)[i],
                         gsl_complex_conjugate( (hptilde->data->data)[i] ) ));
    }
  }
  else if(model_part.compare("hchc") == 0) {
    //fprintf(stdout,"hc conj(hc)\n");
    for (int i=0; i<n; i++) {
      gsl_vector_complex_set(wv, i,
        gsl_complex_mul((hctilde->data->data)[i],
                         gsl_complex_conjugate( (hctilde->data->data)[i] ) ));
    }
  }
  else if(model_part.compare("hpPLUShcSquared") == 0) {
    //fprintf(stdout,"(hc + hp) * conj(hc+hp)\n");
    gsl_vector_complex *wv_conj;
    wv_conj = gsl_vector_complex_alloc(n);
    for(int i=0; i<n; i++) {
      gsl_vector_complex_set(wv, i,
        gsl_complex_add((hctilde->data->data)[i],(hptilde->data->data)[i]));
      gsl_vector_complex_set(wv_conj, i,
        gsl_complex_conjugate( gsl_complex_add((hctilde->data->data)[i],(hptilde->data->data)[i]) ));
    }
    gsl_vector_complex_mul(wv, wv_conj);
    gsl_vector_complex_free(wv_conj);
  }
  else if(model_part.compare("hphc") == 0) {
    //fprintf(stdout,"Im hp hc\n");
    gsl_complex wv_i_real;
    for (int i=0; i<n; i++) {
      const gsl_complex wv_i =
        gsl_complex_mul((hptilde->data->data)[i],
                         gsl_complex_conjugate( (hctilde->data->data)[i] ) ) ;
      double x = GSL_REAL(wv_i);
      GSL_SET_COMPLEX(&wv_i_real, x, 0.0);
      gsl_vector_complex_set(wv, i, wv_i_real);
    }

    // -- HpHc cross term could be very small. Set to zero if too small --//
    gsl_complex size_hphc, size_hp, size_hc;
    gsl_vector_complex *v_tmp;
    v_tmp = gsl_vector_complex_alloc(wv->size);


    // abs_hp, etc. are NOT the norm squared since quadrature weights 
    // are neglected: Thats OK, as we just want a size estimate
    size_hphc = mygsl::EuclideanInner(wv,wv);
    for (int i=0; i<n; i++) {
      gsl_vector_complex_set(v_tmp, i, (hctilde->data->data)[i]);
    }
    size_hc = mygsl::EuclideanInner(v_tmp,v_tmp);
    for (int i=0; i<n; i++) {
      gsl_vector_complex_set(v_tmp, i, (hptilde->data->data)[i]);
    }
    size_hp = mygsl::EuclideanInner(v_tmp,v_tmp);
    const double abs_hp   = gsl_complex_abs(size_hp);
    const double abs_hc   = gsl_complex_abs(size_hc);
    const double abs_hphc = gsl_complex_abs(size_hphc);


    if( (abs_hp + abs_hc) / abs_hphc > 1.e8 ) { // TODO: hardcoded tol is bad 
      fprintf(stdout,"Size hp = %e, size hc = %e, size hphc = %e\n",
              abs_hp, abs_hc, abs_hphc);
      gsl_vector_complex_set_zero(wv);
      size_hphc = mygsl::EuclideanInner(wv,wv);
      fprintf(stdout,"New size hphc = %e\n",gsl_complex_abs(size_hphc));
    }

    gsl_vector_complex_free(v_tmp);

  }
  else {
    std::cerr << "Approximant not supported!" << std::endl;
    exit(1);
  }
}

std::string model_tag2mode_part(const std::string model_tag,
                                const int mode_param_indx,
                                const double *params)
{
  // model_tag must end with "all_parts" or one of the strings
  // specified in lal_waveform_part

  std::string model_part;

  // --- deduce the model_part tag type --- //
  std::size_t found = model_tag.find("all_parts");
  // "all_parts" is a distinct model of higher dimension
  if(found!=std::string::npos) {
    //fprintf(stdout,"all parts model\n");
    model_part = get_waveform_part_tag(params[mode_param_indx]);
  }
  else {
    model_part = get_waveform_part_tag(model_tag);
  }
  //fprintf(stdout,"model part tag %s\n",model_part.c_str());

  // TODO: check that a mode part has been found
  return model_part;

}

std::vector<std::string> get_barycenter_tags(const std::string model_tag){
  std::vector<std::string> x = split_string(model_tag, '_'); // split tag on underscores

  if ( x.size() != 4 && x.size() != 5 ){
    fprintf(stderr, "Model tag should have format \"Barycenter_DET_EPHEM_UNITS\" or \"Barycenter_DET_EPHEM_UNITS_NOSHAPIRO\"\n");
    exit(1);
  }
  // get the last three parts of the vector (i.e. DET, EPHEM, UNITS)
  std::vector<std::string> y(x.begin()+1, x.begin()+4);

  // add whether to include Shapiro delay in model or not
  if ( x.size() == 4 ){ y.push_back("SHAPIRO"); }
  else { y.push_back(x[4]); } 

  // check ephem is DE200, DE405, DE414 or DE421
  if ( strcmp(y[1].c_str(), "DE200") && strcmp(y[1].c_str(), "DE405") &&
       strcmp(y[1].c_str(), "DE414") && strcmp(y[1].c_str(), "DE421") ){
    fprintf(stderr, "Ephemeris must be either \"DE200\", \"DE405\", \"DE414\", or \"DE421\"\n");
    exit(1);
  }

  // check units is either TDB or TCB
  if ( strcmp(y[2].c_str(), "TCB") && strcmp(y[2].c_str(), "TDB" ) ){
    fprintf(stderr, "Time units must be either \"TCB\" or \"TDB\"\n");
    exit(1);
  }

  // check NOSHAPIRO is either SHAPIRO or NOSHAPIRO (i.e. include or remove Shapiro delay from calculations)
  if ( strcmp(y[3].c_str(), "SHAPIRO") && strcmp(y[3].c_str(), "NOSHAPIRO") ){
    fprintf(stderr, "Must specify either \"SHAPIRO\" or \"NOSHAPIRO\"\n");
    exit(1);
  }

  fprintf(stdout, "Detector: \"%s\", ephemeris: \"%s\", time units: \"%s\"\n", y[0].c_str(), y[1].c_str(), y[2].c_str());

  return y;
}

int get_binary_barycenter_tags(const std::string model_tag){
  std::vector<std::string> x = split_string(model_tag, '_'); // split tag on underscores

  if ( x.size() != 1 && x.size() != 2 ){
    fprintf(stderr, "Model tag should have format \"BinaryBarycenter\" or \"BinaryBarycenter_TDOT\"\n");
    exit(1);
  }

  if ( x.size() == 1 ){ return 0; }

  // check for TDOT
  if ( !strcmp(x[1].c_str(), "TDOT") ){ return 1; }
  else {
    fprintf(stderr, "Warning... suffix was not the expected \"TDOT\" value. Defaulting to work with time delays");
    return 0;
  }
}

}
