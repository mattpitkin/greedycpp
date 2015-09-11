// AUTHOR :  Michael Puerrer
//           Michael.Puerrer@ligo.org
//
// DATE: July 1, 2015
//
// PURPOSE: Interface with LAL SEOBNRv2 ROM
//          https://www.lsc-group.phys.uwm.edu/daswg/projects/lalsuite.html
//

#ifndef ROM_SEOBNRv2LAL_HPP
#define ROM_SEOBNRv2LAL_HPP

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


void SEOBNRv2_ROM_SingleSpin_Waveform(gsl_vector_complex *wv,
                    const gsl_vector *fnodes,
                    const double *params)
{

	COMPLEX16FrequencySeries *hptilde = NULL;
	COMPLEX16FrequencySeries *hctilde = NULL;

	int n = fnodes->size;
	double m1SI = params[0] * LAL_MSUN_SI;
	double m2SI = params[1] * LAL_MSUN_SI;
	double chi  = params[2];
	double distance = 100*1e6*LAL_PC_SI;
	double inclination = 0;
	double fRef = 0;
	double phiRef = 0;


	// Copy frequency data into sequence
        const REAL8Sequence *freqs = XLALCreateREAL8Sequence(n);
	for (int i=0; i<n; i++)
    	  freqs->data[i] = gsl_vector_get(fnodes, i);

	/** Compute waveform in LAL format at specified frequencies */
	//int ret = XLALSimIMRSEOBNRv2ROMEqualSpinFrequencySequence(
	int ret = XLALSimIMRSEOBNRv2ROMSingleSpinFrequencySequence(
	  &hptilde,      /**< Output: Frequency-domain waveform h+ */
	  &hctilde,      /**< Output: Frequency-domain waveform hx */
	  freqs,         /**< Frequency points at which to evaluate the waveform (Hz) */
	  phiRef,        /**< Phase at reference time */
	  fRef,          /**< Reference frequency (Hz); 0 defaults to fLow */
	  distance,      /**< Distance of source (m) */
	  inclination,   /**< Inclination of source (rad) */
	  m1SI,          /**< Mass of companion 1 (kg) */
	  m2SI,          /**< Mass of companion 2 (kg) */
	  chi            /**< Dimensionless aligned component spin chi1=chi2 = chi */
	);

	if (ret != XLAL_SUCCESS) {
    	  fprintf(stderr, "Error calling XLALSimIMRSEOBNRv2ROMEqualSpinFrequencySequence().\n");
    	  exit(-1);
  	}
  	XLALDestroyREAL8Sequence((REAL8Sequence *)freqs);


  	// Copy polarization into output buffer
    	for (int i=0; i<n; i++)
	  gsl_vector_complex_set(wv, i, (hptilde->data->data)[i]);

  	XLALDestroyCOMPLEX16FrequencySeries(hptilde);
  	XLALDestroyCOMPLEX16FrequencySeries(hctilde);
}

// TODO: add new DS ROM once its in master
#endif

