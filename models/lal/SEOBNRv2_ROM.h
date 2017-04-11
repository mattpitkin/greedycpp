// AUTHOR :  Michael Puerrer
//           Michael.Puerrer@ligo.org
//
// DATE: 2015-2016
//
// PURPOSE: Interface with LAL SEOBNRv2 ROM and related models
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
	int ret = XLALSimIMRSEOBNRv2ROMEffectiveSpinFrequencySequence(
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

// This routine is interfaced with greedy routine -- returns gsl data type //
void ROM_SEOBNRv2_DS_HI_FullWaveform(
  gsl_vector_complex * &wv, 
  const gsl_vector *fnodes,
  const double *params,
  const std::string model_tag)
{
  // params = {m1,m2,chi1,chi2}
  // parameter list such that (m1(param),m2(param),chi1,chi2) is a unique point in parameter space

  // --- deduce the model_part tag --- //
  std::string model_part =
    lal_help::model_tag2mode_part(model_tag,4,params);


  // Note: We expect masses in units of solar mass 
  // -> use conversion factor 1 in cfg-file!
  double m1SI = params[0] * LAL_MSUN_SI;
  double m2SI = params[1] * LAL_MSUN_SI;
  double chi1 = params[2];
  double chi2 = params[3];
  double distance = 100*1e6*LAL_PC_SI;
  double inclination = 0;
  double fRef = 0;
  double phiRef = 0;


  int n = fnodes->size;
  // Copy frequency data into sequence
  const REAL8Sequence *freqs = XLALCreateREAL8Sequence(n);
  for (int i=0; i<n; i++)
    freqs->data[i] = gsl_vector_get(fnodes, i);

  struct tagCOMPLEX16FrequencySeries *hptilde = NULL;
  struct tagCOMPLEX16FrequencySeries *hctilde = NULL;

  /** Compute waveform in LAL format at specified frequencies */
  int ret = XLALSimIMRSEOBNRv2ROMDoubleSpinHIFrequencySequence(
    &hptilde,      /**< Output: Frequency-domain waveform h+ */
    &hctilde,      /**< Output: Frequency-domain waveform hx */
    freqs,         /**< Frequency points at which to evaluate the waveform (Hz) */
    phiRef,        /**< Phase at reference time */
    fRef,          /**< Reference frequency (Hz); 0 defaults to fLow */
    distance,      /**< Distance of source (m) */
    inclination,   /**< Inclination of source (rad) */
    m1SI,          /**< Mass of companion 1 (kg) */
    m2SI,          /**< Mass of companion 2 (kg) */
    chi1,          /**< Dimensionless aligned component spin 1 */
    chi2,          /**< Dimensionless aligned component spin 2 */
    -1
  );

  if (ret != XLAL_SUCCESS) {
    fprintf(stderr, "Error calling XLALSimIMRSEOBNRv2ROMDoubleSpinHIFrequencySequence().\n");
    exit(-1);
  }

  lal_help::lal_waveform_part(wv,model_part,hptilde,hctilde,n);

  XLALDestroyREAL8Sequence((REAL8Sequence *)freqs);
  XLALDestroyCOMPLEX16FrequencySeries(hptilde);
  XLALDestroyCOMPLEX16FrequencySeries(hctilde);
}


// This routine is interfaced with greedy routine -- returns gsl data type //
void LackeyTidal2013_FullWaveform(
  gsl_vector_complex * &wv, 
  const gsl_vector *fnodes,
  const double *params,
  const std::string model_tag)
{
  // This tidal model is built on top of SEOBNRv2ROM_DS_HI
  // Note that the parameters are different from the plain ROM!

  // params = {mBH,mNS,chiBH,Lambda}
  // parameter list such that (mBH(param),mNS(param),chiBH,Lambda) is a unique point in parameter space

  // --- deduce the model_part tag --- //
  std::string model_part =
    lal_help::model_tag2mode_part(model_tag,4,params);

  // Note: We expect masses in units of solar mass 
  // -> use conversion factor 1 in cfg-file!
  double mBH_SI = params[0] * LAL_MSUN_SI;
  double mNS_SI = params[1] * LAL_MSUN_SI;
  double chiBH = params[2];
  double Lambda = params[3];
  double distance = 100*1e6*LAL_PC_SI;
  double inclination = 0;
  double fRef = 0;
  double phiRef = 0;


  int n = fnodes->size;
  // Copy frequency data into sequence
  const REAL8Sequence *freqs = XLALCreateREAL8Sequence(n);
  for (int i=0; i<n; i++)
    freqs->data[i] = gsl_vector_get(fnodes, i);

  struct tagCOMPLEX16FrequencySeries *hptilde = NULL;
  struct tagCOMPLEX16FrequencySeries *hctilde = NULL;

  // Compute waveform in LAL format at specified frequencies //
  int ret = XLALSimIMRLackeyTidal2013FrequencySequence(
    &hptilde,      //< Output: Frequency-domain waveform h+ //
    &hctilde,      //< Output: Frequency-domain waveform hx //
    freqs,         //< Frequency points at which to evaluate the waveform (Hz) //
    phiRef,        //< Phase at reference time //
    fRef,          //< Reference frequency (Hz); 0 defaults to fLow //
    distance,      //< Distance of source (m) //
    inclination,   //< Inclination of source (rad) //
    mBH_SI,        //< Mass of black hole (kg) //
    mNS_SI,        //< Mass of neutron star (kg) //
    chiBH,         //< Dimensionless aligned component spin of the BH //
    Lambda         //< Dimensionless tidal deformability (Eq 1  of Lackey et al) //
  );

  if (ret != XLAL_SUCCESS) {
    fprintf(stderr, "Error calling XLALSimIMRLackeyTidal2013FrequencySequence().\n");
    exit(-1);
  }

  lal_help::lal_waveform_part(wv,model_part,hptilde,hctilde,n);


  XLALDestroyREAL8Sequence((REAL8Sequence *)freqs);
  XLALDestroyCOMPLEX16FrequencySeries(hptilde);
  XLALDestroyCOMPLEX16FrequencySeries(hctilde);
}

#endif

