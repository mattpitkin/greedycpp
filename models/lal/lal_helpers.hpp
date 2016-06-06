// AUTHOR :  Scott Field 
//           sfield@astro.cornell.edu
//
// DATE: May 30, 2016
//
// PURPOSE: various lal helper routines

#ifndef lal_helpers_h
#define lal_helpers_h

#include <string>

#include <lal/FrequencySeries.h>
#include "../../code/gsl_helper_functions.hpp"

namespace lal_help {

// deduce waveform part from full model tag's substring
// Requires model tag to be of the form NAME_XXX
// valid tags are those cases listed in get_waveform_part_tag
std::string get_waveform_part_tag(const std::string model_tag);

// an "all_parts" model adds an extra parameter to denote the model
// part. convert this value, model_type, into a model_part string id
std::string get_waveform_part_tag(const double model_type);

// Copy polarization into output buffer
void lal_waveform_part(gsl_vector_complex *wv,
                       const std::string model_part,
                       const COMPLEX16FrequencySeries *hptilde,
                       const COMPLEX16FrequencySeries *hctilde,
                       const int n);

// Decide if mode part in models tag or "all_parts" model and
// return model_part string id
std::string model_tag2mode_part(const std::string model_tag,
                                const int mode_param_indx,
                                const double *params);

};

#endif
