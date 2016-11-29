extern "C"{
#ifdef __cplusplus
  // Workaround to allow C++ programs to use stdint.h macros specified in the C99 standard that aren't in the C++ standard.
  #define __STDC_CONSTANT_MACROS
  #ifdef _STDINT_H
    #undef _STDINT_H
  #endif
  #include <stdint.h>
#endif
}

//#include <iostream>
#include <stdio.h>
#include <math.h>

#include <lal/LALConstants.h>
#include <lal/LALSimIMR.h>
//#include <stdint.h>
#include <lal/Sequence.h>
#include <lal/FrequencySeries.h>

#include "lal_helpers.hpp"

/*----------------------------------------------------------------------------
	Calls a modified version of XLALSimIMRPhenomP that allows the waveform
	to be computed for an array of arbitrary ordered frequencies.
 -----------------------------------------------------------------------------*/


void PhenP_Waveform(gsl_vector_complex *wv,
                    const gsl_vector *fnodes,
                    const double *params,
                    const std::string model_tag)
{


  // --- deduce the model_part tag --- //
  std::string model_part =
    lal_help::model_tag2mode_part(model_tag,7,params);


  COMPLEX16FrequencySeries *hptilde = NULL;
  COMPLEX16FrequencySeries *hctilde = NULL;

  int n = fnodes->size;
  const REAL8 m1_Msun = params[0];
  const REAL8 m2_Msun = params[1];
  const REAL8 m1_SI = m1_Msun * LAL_MSUN_SI;
  const REAL8 m2_SI = m2_Msun * LAL_MSUN_SI;
  const REAL8 chi1L = params[2];
  const REAL8 chi2L = params[3];
  const REAL8 chip = params[4];
  const REAL8 thetaJ = params[5];
  const REAL8 distance = 1;
  const REAL8 phic = 0;
  const REAL8 f_ref = 40;
  const REAL8 alpha0 = params[6];
  IMRPhenomP_version_type version_flag = IMRPhenomPv2_V;

  // use XLALSimIMRPhenomPFrequencySequence with frequency sequence //
  const REAL8Sequence *freqs = XLALCreateREAL8Sequence(n);
  for (int i=0; i<n; i++) {
    freqs->data[i] = gsl_vector_get(fnodes, i);
  }

  int ret = XLALSimIMRPhenomPFrequencySequence(
    &hptilde,   //< Output: Frequency-domain waveform h+ //
    &hctilde,   //< Output: Frequency-domain waveform hx //
    freqs,           //< Frequency points at which to evaluate the waveform (Hz) //
    chi1L,                  //< Effective aligned spin //
    chi2L,
    chip,                     //< Effective spin in the orbital plane //
    thetaJ,
    m1_SI,                      //< Symmetric mass-ratio //
    m2_SI,                   //< Angle between J0 and line of sight (z-direction) //
    distance,                  //< Total mass of binary (kg) //
    alpha0,                 //< Distance of source (m) //
    phic,                   //< Initial value of alpha angle (azimuthal precession angle) //
    f_ref,                     //< Orbital phase at the peak of the underlying non precessing model (rad) //
    version_flag,                    //< Reference frequency //
    NULL);

  if (ret != XLAL_SUCCESS) {
    fprintf(stderr, "Error calling XLALSimIMRPhenomPFrequencySequence().\n");
    exit(-1);
  }


  lal_help::lal_waveform_part(wv,model_part,hptilde,hctilde,n);

  XLALDestroyREAL8Sequence((REAL8Sequence *)freqs);
  XLALDestroyCOMPLEX16FrequencySeries(hptilde);
  XLALDestroyCOMPLEX16FrequencySeries(hctilde);
}



// phenomP version 1 with older interface 
// use LAL git hash a27aef328a77a5de5434c27d22f812bfe369c7a8
/*void PhenP_Waveform(gsl_vector_complex *wv,
                    const gsl_vector *fnodes,
                    const double *params,
                    const std::string model_tag)
{


  // --- deduce the model_part tag --- //
  std::string model_part;
  // "all_parts" is a distinct model of higher dimension (i.e. the
  // 7th describes the "part") and it needs its own model tag
  if(model_tag.compare("PhenomP_all_parts") == 0) {
    fprintf(stdout,"all parts model\n");
    // 7th parameter params[6] 0 or 1 for plus or cross mode, 2-4 for quad part
    model_part = get_waveform_part_tag(params[6]);
  }
  else {
    model_part = get_waveform_part_tag(model_tag);
  }
  fprintf(stdout,"model part tag %s\n",model_part.c_str());



  COMPLEX16FrequencySeries *hptilde = NULL;
  COMPLEX16FrequencySeries *hctilde = NULL;

  int n = fnodes->size;
  const REAL8 Mtot_Msun = params[0];
  const REAL8 eta = params[1];
  const REAL8 Mtot_SI = Mtot_Msun * LAL_MSUN_SI; 
  const REAL8 chi_eff = params[2];
  const REAL8 chip = params[3];
  const REAL8 thetaJ = params[4];
  const REAL8 distance = 1;
  const REAL8 phic = 0;
  const REAL8 f_ref = 40;
  const REAL8 alpha0 = params[5];


  // use XLALSimIMRPhenomPFrequencySequence with frequency sequence //
  const REAL8Sequence *freqs = XLALCreateREAL8Sequence(n);
  for (int i=0; i<n; i++) {
    freqs->data[i] = gsl_vector_get(fnodes, i);
  }

  int ret = XLALSimIMRPhenomPFrequencySequence(
    &hptilde,   //< Output: Frequency-domain waveform h+ //
    &hctilde,   //< Output: Frequency-domain waveform hx //
    freqs,           //< Frequency points at which to evaluate the waveform (Hz) //
    chi_eff,                  //< Effective aligned spin //
    chip,                     //< Effective spin in the orbital plane //
    eta,                      //< Symmetric mass-ratio //
    thetaJ,                   //< Angle between J0 and line of sight (z-direction) //
    Mtot_SI,                  //< Total mass of binary (kg) //
    distance,                 //< Distance of source (m) //
    alpha0,                   //< Initial value of alpha angle (azimuthal precession angle) //
    phic,                     //< Orbital phase at the peak of the underlying non precessing model (rad) //
    f_ref);                    //< Reference frequency //

  if (ret != XLAL_SUCCESS) {
    fprintf(stderr, "Error calling XLALSimIMRPhenomPFrequencySequence().\n");
    exit(-1);
  }
  XLALDestroyREAL8Sequence((REAL8Sequence *)freqs);

  lal_waveform_part(wv,model_part,hptilde,hctilde,n);

  XLALDestroyCOMPLEX16FrequencySeries(hptilde);
  XLALDestroyCOMPLEX16FrequencySeries(hctilde);
}*/


/*
// KEEP -- old interface which agrees with new frequency series //
void PhenP_Waveform(gsl_vector_complex *wv, const gsl_vector *fnodes, double *params, const char *plus_cross_flag)
{

	COMPLEX16FrequencySeries *hptilde = NULL;
	COMPLEX16FrequencySeries *hctilde = NULL;

	int n = fnodes->size;
	int offset = 0; // XLALSimIMRPhenomP pads f<fmin with zeros

	const REAL8 Mtot_Msun = params[0];
	const REAL8 eta = params[1];
	const REAL8 Mtot_SI = Mtot_Msun * LAL_MSUN_SI; 
	const REAL8 chi_eff = params[2];
	const REAL8 chip = params[3];
	const REAL8 thetaJ = params[4];
	const REAL8 distance = 1;
	const REAL8 phic = 0;
	const REAL8 f_ref = 40;
	const REAL8 alpha0 = params[5];


	// used XLALSimIMRPhenomP with fixed deltaF //
	const REAL8 f_min  = gsl_vector_get(fnodes,0);
	const REAL8 deltaF = gsl_vector_get(fnodes,1) - gsl_vector_get(fnodes,0);
	const REAL8 f_max  = gsl_vector_get(fnodes,fnodes->size-1) + deltaF;
	double deltaF_double =  gsl_vector_get(fnodes,1) - gsl_vector_get(fnodes,0);
	double fmin_double   = gsl_vector_get(fnodes,0);

	XLALSimIMRPhenomP(
	  &hptilde,   //< Output: Frequency-domain waveform h+ /
	  &hctilde,   //< Output: Frequency-domain waveform hx //
	  chi_eff,                  //< Effective aligned spin //
	  chip,                     //< Effective spin in the orbital plane //
	  eta,                      //< Symmetric mass-ratio //
	  thetaJ,                   //< Angle between J0 and line of sight (z-direction) //
	  Mtot_SI,                  //< Total mass of binary (kg) //
	  distance,                 //< Distance of source (m) //
	  alpha0,                   //< Initial value of alpha angle //
	  phic,                     //< Orbital phase at the peak of the underlying non precessing model (rad) //
	  deltaF,                   //< Sampling frequency (Hz) //
	  f_min,                    //< Starting GW frequency (Hz) //
	  f_max,                    //< End frequency; 0 defaults to ringdown cutoff freq //
	  f_ref);                    //< Reference frequency //

	// check this is indeed the offset -- data[offset] == 0, data[offset+1]!=0 //
	// LAL code pads frequencies 0 through fmin-deltaF with zero data //
	offset = (int) ((int) (fmin_double/deltaF_double));
	if( creal((hptilde->data->data)[offset-1]) !=0 || cimag((hptilde->data->data)[offset-1]) !=0) {
	  fprintf(stderr, "nonzero data.\n");
    	  exit(-1);
	}
	if( creal((hptilde->data->data)[offset]) ==0 || cimag((hptilde->data->data)[offset]) ==0) {
	  fprintf(stderr, "zero data.\n");
    	  exit(-1);
	}
	offset = offset; // offset used below whe copying into output buffer 

  	// Copy polarization into output buffer
	// used for both XLALSimIMRPhenomPFrequencySequence and XLALSimIMRPhenomP
	if(strcmp(plus_cross_flag,"PhenomP_plus") == 0) {
    	  for (int i=offset; i<(offset+n); i++)
	    gsl_vector_complex_set(wv, i-offset, (hptilde->data->data)[i]);
	}
  	else if(strcmp(plus_cross_flag,"PhenomP_cross") == 0) {
    	  for (int i=offset; i<(offset+n); i++)
      	    gsl_vector_complex_set(wv, i-offset, (hctilde->data->data)[i]);
  	}

  	XLALDestroyCOMPLEX16FrequencySeries(hptilde);
  	XLALDestroyCOMPLEX16FrequencySeries(hctilde);
}*/
