//#include <iostream>
#include <stdio.h>
#include <math.h>
#include <lal/LALConstants.h>
#include <lal/LALSimIMR.h>
#include <stdint.h>
#include <lal/FrequencySeries.h>
/*----------------------------------------------------------------------------
	Calls a modified version of XLALSimIMRPhenomP that allows the waveform
	to be computed for an array of arbitrary ordered frequencies.
 -----------------------------------------------------------------------------*/


void PhenP_Waveform(gsl_vector_complex *wv, const gsl_vector *fnodes, double *params, const char *plus_cross_flag)
{

	COMPLEX16FrequencySeries *hptilde = NULL;
	COMPLEX16FrequencySeries *hctilde = NULL;


	const REAL8 Mtot_Msun = params[0];
	const REAL8 eta = params[1];
	const REAL8 Mtot_SI = Mtot_Msun * LAL_MSUN_SI; 
	const REAL8 chi_eff = params[2];
	const REAL8 chip = params[3];
	const REAL8 thetaJ = params[4];
	const REAL8 phiJ = 0;
	const REAL8 distance = 1;
	const REAL8 phic = 0;
	const REAL8 f_ref = 40;
	const REAL8 alpha0 = params[5];

	const REAL8 f_min  = gsl_vector_get(fnodes,0);
	const REAL8 f_max  = gsl_vector_get(fnodes,fnodes->size);
	const REAL8 deltaF = gsl_vector_get(fnodes,1) - gsl_vector_get(fnodes,0);

	XLALSimIMRPhenomP(
	  &hptilde,   /**< Output: Frequency-domain waveform h+ */
	  &hctilde,   /**< Output: Frequency-domain waveform hx */
	  chi_eff,                  /**< Effective aligned spin */
	  chip,                     /**< Effective spin in the orbital plane */
	  eta,                      /**< Symmetric mass-ratio */
	  thetaJ,                   /**< Angle between J0 and line of sight (z-direction) */
	  phiJ,                     /**< Angle of J0 in the plane of the sky */
	  Mtot_SI,                  /**< Total mass of binary (kg) */
	  distance,                 /**< Distance of source (m) */
	  alpha0,                   /**< Initial value of alpha angle */
	  phic,                     /**< Orbital phase at the peak of the underlying non precessing model (rad) */
	  deltaF,                   /**< Sampling frequency (Hz) */
	  f_min,                    /**< Starting GW frequency (Hz) */
	  f_max,                    /**< End frequency; 0 defaults to ringdown cutoff freq */
	  f_ref);                    /**< Reference frequency */
	/* removed Rory's XLALSimIMRPhenomP_forRB which takes in arbitrary frequencies (future lal patch)
        XLALSimIMRPhenomP_forRB(
          &hptilde,   //< Output: Frequency-domain waveform h+ /
          &hctilde,   //< Output: Frequency-domain waveform hx 
          chi_eff,                  //< Effective aligned spin 
          chip,                     //< Effective spin in the orbital plane 
          eta,                      /< Symmetric mass-ratio 
          thetaJ,                   /< Angle between J0 and line of sight (z-direction) /
          phiJ,                     /< Angle of J0 in the plane of the sky /
          Mtot_SI,                  /< Total mass of binary (kg) /
          distance,                 /< Distance of source (m) /
          alpha0,                   /< Initial value of alpha angle /
          phic,                     /< Orbital phase at the peak of the underlying non precessing model (rad) /
          f_ref,
          fnodes);*/


	gsl_complex zM;

	if(strcmp(plus_cross_flag,"PhenomP_plus") == 0){

  		for(int cols = 0; cols < fnodes->size; cols++)
   		{
        		GSL_SET_COMPLEX(&zM, creal(hptilde->data->data[cols]), cimag(hptilde->data->data[cols]));
        		gsl_vector_complex_set(wv,cols,zM);
    		}
	}

	else if(strcmp(plus_cross_flag,"PhenomP_cross") == 0){

  		for(int cols = 0; cols < fnodes->size; cols++)
   		{
        		GSL_SET_COMPLEX(&zM, creal(hctilde->data->data[cols]), cimag(hctilde->data->data[cols]));
        		gsl_vector_complex_set(wv,cols,zM);
    		}
	}

	XLALDestroyCOMPLEX16FrequencySeries(hptilde);
	XLALDestroyCOMPLEX16FrequencySeries(hctilde);
}


void hp_hc_hphc(gsl_vector_complex *wv, const gsl_vector *fnodes, double *params, const char *plus_cross_flag)
{

	COMPLEX16FrequencySeries *hptilde = NULL;
	COMPLEX16FrequencySeries *hctilde = NULL;

	const REAL8 Mtot_Msun = params[0];
	const REAL8 eta = params[1];
	const REAL8 Mtot_SI = Mtot_Msun * LAL_MSUN_SI; 
	const REAL8 chi_eff = params[2];
	const REAL8 chip = params[3];
	const REAL8 thetaJ = params[4];
	const REAL8 phiJ = 0;
	const REAL8 distance = 1;
	const REAL8 phic = 0;
	const REAL8 f_ref = 40;
	const REAL8 alpha0 = params[5];

	const REAL8 f_min  = gsl_vector_get(fnodes,0);

	const REAL8 f_max  = gsl_vector_get(fnodes,fnodes->size-1);
	const REAL8 deltaF = gsl_vector_get(fnodes,1) - gsl_vector_get(fnodes,0);

	XLALSimIMRPhenomP(
	  &hptilde,   /**< Output: Frequency-domain waveform h+ */
	  &hctilde,   /**< Output: Frequency-domain waveform hx */
	  chi_eff,                  /**< Effective aligned spin */
	  chip,                     /**< Effective spin in the orbital plane */
	  eta,                      /**< Symmetric mass-ratio */
	  thetaJ,                   /**< Angle between J0 and line of sight (z-direction) */
	  phiJ,                     /**< Angle of J0 in the plane of the sky */
	  Mtot_SI,                  /**< Total mass of binary (kg) */
	  distance,                 /**< Distance of source (m) */
	  alpha0,                   /**< Initial value of alpha angle */
	  phic,                     /**< Orbital phase at the peak of the underlying non precessing model (rad) */
	  deltaF,                   /**< Sampling frequency (Hz) */
	  f_min,                    /**< Starting GW frequency (Hz) */
	  f_max,                    /**< End frequency; 0 defaults to ringdown cutoff freq */
	  f_ref);                    /**< Reference frequency */
	/* removed Rory's XLALSimIMRPhenomP_forRB which takes in arbitrary frequencies (future lal patch)
        XLALSimIMRPhenomP_forRB(
          &hptilde,   //< Output: Frequency-domain waveform h+ /
          &hctilde,   //< Output: Frequency-domain waveform hx 
          chi_eff,                  //< Effective aligned spin 
          chip,                     //< Effective spin in the orbital plane 
          eta,                      /< Symmetric mass-ratio 
          thetaJ,                   /< Angle between J0 and line of sight (z-direction) /
          phiJ,                     /< Angle of J0 in the plane of the sky /
          Mtot_SI,                  /< Total mass of binary (kg) /
          distance,                 /< Distance of source (m) /
          alpha0,                   /< Initial value of alpha angle /
          phic,                     /< Orbital phase at the peak of the underlying non precessing model (rad) /
          f_ref,
          fnodes);*/

	gsl_complex zM;

	if(strcmp(plus_cross_flag,"hphp") == 0){
  		for(int cols = 0; cols < fnodes->size; cols++)
   		{
                      //std::cout << cols << std::endl;
                      //std::cout << "real: " << creal(hptilde->data->data[cols]) << "imag: " << cimag(hptilde->data->data[cols]) << std::endl;
        		GSL_SET_COMPLEX(&zM, creal(hptilde->data->data[cols])*creal(hptilde->data->data[cols]) + cimag(hptilde->data->data[cols])*cimag(hptilde->data->data[cols]), 0);
        		gsl_vector_complex_set(wv,cols,zM);
    		}
	}

	else if(strcmp(plus_cross_flag,"hchc") == 0){

  		for(int cols = 0; cols < fnodes->size; cols++)
   		{
        		GSL_SET_COMPLEX(&zM, creal(hctilde->data->data[cols])*creal(hctilde->data->data[cols]) + cimag(hctilde->data->data[cols])*cimag(hctilde->data->data[cols]), 0.);
        		gsl_vector_complex_set(wv,cols,zM);
    		}
	}

	else if(strcmp(plus_cross_flag,"hphc") == 0){

                for(int cols = 0; cols < fnodes->size; cols++)
                {
                        GSL_SET_COMPLEX(&zM, creal(hptilde->data->data[cols])*creal(hctilde->data->data[cols]) + cimag(hptilde->data->data[cols])*cimag(hctilde->data->data[cols]), 0.);
                        gsl_vector_complex_set(wv,cols,zM);
                }
        }

	XLALDestroyCOMPLEX16FrequencySeries(hptilde);
	XLALDestroyCOMPLEX16FrequencySeries(hctilde);
}
