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


void PhenP_Waveform(gsl_vector_complex *wv, const gsl_vector *fnodes, const double *params, const char *plus_cross_flag)
{

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
	for (int i=0; i<n; i++)
    	  freqs->data[i] = gsl_vector_get(fnodes, i);

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


  	// Copy polarization into output buffer
	if(strcmp(plus_cross_flag,"PhenomP_plus") == 0) {
    	  for (int i=0; i<n; i++)
	    gsl_vector_complex_set(wv, i, (hptilde->data->data)[i]);
	}
  	else if(strcmp(plus_cross_flag,"PhenomP_cross") == 0) {
    	  for (int i=0; i<n; i++)
      	    gsl_vector_complex_set(wv, i, (hctilde->data->data)[i]);
  	}
  	else if(strcmp(plus_cross_flag,"PhenomP_hphp") == 0) {
    	  for (int i=0; i<n; i++)
      	    gsl_vector_complex_set(wv, i,
              gsl_complex_mul((hptilde->data->data)[i], gsl_complex_conjugate( (hptilde->data->data)[i] ) ));
  	}
  	else if(strcmp(plus_cross_flag,"PhenomP_hchc") == 0) {
    	  for (int i=0; i<n; i++)
      	    gsl_vector_complex_set(wv, i,
              gsl_complex_mul((hctilde->data->data)[i], gsl_complex_conjugate( (hctilde->data->data)[i] ) ));
  	}
  	else if(strcmp(plus_cross_flag,"PhenomP_hphc") == 0) {
    	  std::cerr << "Double check notes for correct expression" << std::endl;
	  exit(1);
	  //for (int i=0; i<n; i++)
      	  //  gsl_vector_complex_set(wv, i, gsl_complex_mul((hptilde->data->data)[i], (hctilde->data->data)[i]));
  	}
	else {
	  std::cerr << "Approximant not supported!" << std::endl;
	  exit(1);
	}

  	// Copy products of polarizations into output buffer
	/*gsl_complex zM;

	else if(strcmp(plus_cross_flag,"PhenomP_hphc") == 0){

                for(int cols = 0; cols < fnodes->size; cols++)
                {
                        GSL_SET_COMPLEX(&zM, creal(hptilde->data->data[cols])*creal(hctilde->data->data[cols]) + cimag(hptilde->data->data[cols])*cimag(hctilde->data->data[cols]), 0.);
                        gsl_vector_complex_set(wv,cols,zM);
                }
        }*/


  	XLALDestroyCOMPLEX16FrequencySeries(hptilde);
  	XLALDestroyCOMPLEX16FrequencySeries(hctilde);
}


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
