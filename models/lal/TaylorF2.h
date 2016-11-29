// AUTHOR :  Michael Puerrer
//           Michael.Puerrer@ligo.org
//
// DATE: July 1, 2015
//
// PURPOSE: Interface with LAL TaylorF2
//          https://www.lsc-group.phys.uwm.edu/daswg/projects/lalsuite.html
//

#ifndef TaylorF2_LAL_HPP
#define TaylorF2_LAL_HPP

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

#include <stdio.h>
#include <math.h>
#include <lal/LALConstants.h>
#include <lal/LALSimIMR.h>
#include <lal/Sequence.h>
#include <lal/FrequencySeries.h>


void TaylorF2_LAL_Waveform(gsl_vector_complex *wv,
                    const gsl_vector *fnodes,
                    const double *params)
{

	COMPLEX16FrequencySeries *hptilde = NULL;

	int n = fnodes->size;
	double m1SI = params[0] * LAL_MSUN_SI;
	double m2SI = params[1] * LAL_MSUN_SI;
	double chi1 = 0; // Don't include spins as parameters for now
	double chi2 = 0;
	double distance = 100*1e6*LAL_PC_SI;
	double inclination = 0;
	double fRef = 0;
	double phiRef = 0;
	int phaseO = -1;
	int amplitudeO = -1;
	double quadparam1 = 1.0; // quadrupole deformation parameter of body 1 (dimensionless, 1 for BH)
	double quadparam2 = 1.0; // quadrupole deformation parameter of body 1 (dimensionless, 1 for BH)
	double lambda1 = 0.0; // (tidal deformation of body 1)/(mass of body 1)^5
	double lambda2 = 0.0; // (tidal deformation of body 1)/(mass of body 1)^5

	// Copy frequency data into sequence
        const REAL8Sequence *freqs = XLALCreateREAL8Sequence(n);
	for (int i=0; i<n; i++)
    	  freqs->data[i] = gsl_vector_get(fnodes, i);

	/** Compute waveform in LAL format at specified frequencies */
	int ret = XLALSimInspiralTaylorF2Core(
	  &hptilde, 
	  freqs, 
	  phiRef,
	  m1SI, m2SI, 
	  chi1, chi2, 
	  fRef,
	  inclination, 
	  distance,
	  quadparam1,
	  quadparam2,
	  lambda1,
	  lambda2,
	  LAL_SIM_INSPIRAL_SPIN_ORDER_ALL,
	  LAL_SIM_INSPIRAL_TIDAL_ORDER_0PN,
          phaseO, amplitudeO,NULL);

	if (ret != XLAL_SUCCESS) {
    	  fprintf(stderr, "Error calling XLALSimInspiralTaylorF2Core().\n");
    	  exit(-1);
  	}
  	XLALDestroyREAL8Sequence((REAL8Sequence *)freqs);


  	// Copy polarization into output buffer
    	for (int i=0; i<n; i++)
	  gsl_vector_complex_set(wv, i, (hptilde->data->data)[i]);

  	XLALDestroyCOMPLEX16FrequencySeries(hptilde);
}

#endif

