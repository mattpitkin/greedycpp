// Author: Matt Pitkin (2016)
// Interface with lalpulsar solar system barycentering codes

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
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/SFTutils.h>

// global variables for barycentering (so ephemeris files are only read in once)
EphemerisData *edat = NULL;
TimeCorrectionData *tdat = NULL;

// set Earth and Sun solar system files (hardcode to DE405)
char earthfile[256] = "earth00-19-DE405.dat.gz";
char sunfile[256] = "sun00-19-DE405.dat.gz";
char timecorrfile[256] = "te405_2000-2019.dat.gz"; // TCB time correction file

void Barycenter_Waveform(gsl_vector_complex *wv,
                    const gsl_vector *timestamps,
                    const double *params,
                    const std::string model_tag){
  // deduce the detector
  LALDetector det;
  if ( strstr(model_tag.c_str(), "H1") != NULL ){ 
    det = *XLALGetSiteInfo( "H1" );
  }
  else if ( strstr(model_tag.c_str(), "L1") != NULL ){ 
    det = *XLALGetSiteInfo( "L1" );
  }
  else if ( strstr(model_tag.c_str(), "V1") != NULL ){ 
    det = *XLALGetSiteInfo( "V1" );
  }
  else{
    fprintf(stderr, "Error... detector not found.\n");
    exit(-1);
  }

  int n = timestamps->size;
  REAL8 ra = params[0];  // right ascension
  REAL8 dec = params[1]; // declination
  
  // variables for calculating barycenter time delay
  BarycenterInput baryinput;
  EarthState earth;
  EmissionTime emit;
  
  if ( !edat ){ edat = XLALInitBarycenter( earthfile, sunfile ); }
  if ( !tdat ){ tdat = XLALInitTimeCorrections( timecorrfile ); }

  /* set up location of detector */
  baryinput.site.location[0] = det.location[0]/LAL_C_SI;
  baryinput.site.location[1] = det.location[1]/LAL_C_SI;
  baryinput.site.location[2] = det.location[2]/LAL_C_SI;
  
  // set source position
  baryinput.dInv = 0.;
  baryinput.delta = dec;
  baryinput.alpha = ra;
  
  for ( int i=0; i < n; i++ ){
    XLALGPSSetREAL8(&baryinput.tgps, gsl_vector_get(timestamps, i));
    XLALBarycenterEarthNew( &earth, &baryinput.tgps, edat, tdat, TIMECORRECTION_TCB );
    XLALBarycenter( &emit, &baryinput, &earth );
    gsl_complex emitdt;
    GSL_SET_COMPLEX(&emitdt, emit.deltaT, 0.);

    // fill in the output training buffer
    gsl_vector_complex_set(wv, i, emitdt);
  }
}