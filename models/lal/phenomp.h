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
/*----------------------------------------------------------------------------
	Calls a modified version of XLALSimIMRPhenomP that allows the waveform
	to be computed for an array of arbitrary ordered frequencies.
 -----------------------------------------------------------------------------*/


// deduce waveform part from full model tag's substring
// Requires model tag to be of the form NAME_XXX
// where XXX is plus, cross, hphp, hchc or hphc
std::string get_waveform_part_tag(const std::string model_tag)
{
  std::size_t pos  = model_tag.find("_");
  std::string part = model_tag.substr(pos+1);
  return part;
}

// an "all_parts" model adds an extra parameter to denote the model
// part. convert this value, model_type, into a model_part string id
std::string get_waveform_part_tag(const double model_type)
{
  std::string model_part;
  if( std::abs( model_type - 0) < 1.e-10) {
    fprintf(stdout,"Plus with 7th param %f\n",model_type);
    //PhenP_Waveform(wv,fnodes,params,"PhenomP_plus");
    model_part = "plus";
  }
  else if( std::abs( model_type - 1) < 1.e-10) {
    fprintf(stdout,"cross with 7th param %f\n",model_type);
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
  else {
    std::cerr << "PhenomP all parts -- unknown part" << std::endl;
    exit(1);
  }
  return model_part;
}

// Copy polarization into output buffer
void lal_waveform_part(gsl_vector_complex *wv,
                       const std::string model_part,
                       const COMPLEX16FrequencySeries *hptilde,
                       const COMPLEX16FrequencySeries *hctilde,
                       const int n)
{

  if(model_part.compare("plus") == 0) {
    fprintf(stdout,"Im plus\n");
    for (int i=0; i<n; i++) {
      gsl_vector_complex_set(wv, i, (hptilde->data->data)[i]);
    }
  }
  else if(model_part.compare("cross") == 0) {
    fprintf(stdout,"Im cross\n");
    for (int i=0; i<n; i++) {
      gsl_vector_complex_set(wv, i, (hctilde->data->data)[i]);
    }
  }
  else if(model_part.compare("hphp") == 0) {
    fprintf(stdout,"Im hp hp\n");
    for (int i=0; i<n; i++) {
      gsl_vector_complex_set(wv, i,
        gsl_complex_mul((hptilde->data->data)[i],
                         gsl_complex_conjugate( (hptilde->data->data)[i] ) ));
    }
  }
  else if(model_part.compare("hchc") == 0) {
    fprintf(stdout,"Im hc hc\n");
    for (int i=0; i<n; i++) {
      gsl_vector_complex_set(wv, i,
        gsl_complex_mul((hctilde->data->data)[i],
                         gsl_complex_conjugate( (hctilde->data->data)[i] ) ));
    }
  }
  else if(model_part.compare("hphc") == 0) {
    fprintf(stdout,"Im hp hc\n");
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



void PhenP_Waveform(gsl_vector_complex *wv,
                    const gsl_vector *fnodes,
                    const double *params,
                    const std::string model_tag)
{


  // --- deduce the model_part tag --- //
  std::string model_part;
  // "all_parts" is a distinct model of higher dimension
  if(model_tag.compare("PhenomP_all_parts") == 0) {
    fprintf(stdout,"all parts model\n");
    model_part = get_waveform_part_tag(params[7]);
  }
  else {
    model_part = get_waveform_part_tag(model_tag);
  }
  fprintf(stdout,"model part tag %s\n",model_part.c_str());


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
    version_flag);                    //< Reference frequency //

  if (ret != XLAL_SUCCESS) {
    fprintf(stderr, "Error calling XLALSimIMRPhenomPFrequencySequence().\n");
    exit(-1);
  }
  XLALDestroyREAL8Sequence((REAL8Sequence *)freqs);

  lal_waveform_part(wv,model_part,hptilde,hctilde,n);

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
