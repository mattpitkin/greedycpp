// Author: Matt Pitkin (2016)
// Interface with lalpulsar solar system and binary system barycentering codes

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
#include <lal/BinaryPulsarTiming.h>

#include <omp.h>

// global variables for barycentering (so ephemeris files are only read in once)
EphemerisData *edat = NULL;
TimeCorrectionData *tdat = NULL;
LALDetector det;
TimeCorrectionType ttype;
int noshapiro = 0;

void Barycenter_Waveform(gsl_vector_complex *wv,
                    const gsl_vector *timestamps,
                    const double *params,
                    const std::string model_tag){
  int n = timestamps->size;
  REAL8 ra = params[0];  // right ascension
  REAL8 dec = params[1]; // declination

  // variables for calculating barycenter time delay
  EarthState earth;
  EmissionTime emit;
  BarycenterInput baryinput;

  if ( !edat && !tdat ){
    #ifdef USE_OPENMP
    #pragma omp master
    {
    #endif
    // deduce the detector, ephemeris and time units
    std::vector<std::string> vals = lal_help::get_barycenter_tags(model_tag);

    // TODO: add more detectors (and radio telescopes) in the future
    if ( strstr(vals[0].c_str(), "H1") != NULL ){ 
      det = *XLALGetSiteInfo( "H1" );
    }
    else if ( strstr(vals[0].c_str(), "L1") != NULL ){ 
      det = *XLALGetSiteInfo( "L1" );
    }
    else if ( strstr(vals[0].c_str(), "V1") != NULL ){ 
      det = *XLALGetSiteInfo( "V1" );
    }
    else{
      fprintf(stderr, "Error... detector not found.\n");
      exit(-1);
    }

    // set earth and sun ephemeris files
    char earthfile[256], sunfile[256];
    snprintf(earthfile, sizeof(char)*256, "earth00-19-%s.dat.gz", vals[1].c_str());
    snprintf(sunfile, sizeof(char)*256, "sun00-19-%s.dat.gz", vals[1].c_str());
    edat = XLALInitBarycenter( earthfile, sunfile );

    // set time units file
    char timecorrfile[256];
    if ( !strcmp("TCB", vals[2].c_str()) ){
      snprintf(timecorrfile, sizeof(char)*256, "te405_2000-2019.dat.gz");
      ttype = TIMECORRECTION_TCB;
    }
    else if ( !strcmp("TDB", vals[2].c_str()) ){
      snprintf(timecorrfile, sizeof(char)*256, "tdb_2000-2019.dat.gz");
      ttype = TIMECORRECTION_TDB;
    }
    else{
      fprintf(stderr, "Error... time system must be \"TDB\" or \"TCB\".\n");
      exit(-1);
    }
    tdat = XLALInitTimeCorrections( timecorrfile );
    
    if ( !strcmp("NOSHAPIRO", vals[3].c_str()) ){
      noshapiro = 1; // do not include Shapiro delay in barycenter calculation
    }
    #ifdef USE_OPENMP
    }
    #pragma omp barrier
    #endif
  }

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
    XLALBarycenterEarthNew( &earth, &baryinput.tgps, edat, tdat, ttype );
    XLALBarycenter( &emit, &baryinput, &earth );
    gsl_complex emitdt;
    if ( noshapiro ){ // remove Shapiro delay
      GSL_SET_COMPLEX(&emitdt, emit.deltaT+emit.shapiro, 0.); // add Shapiro delay as it is subtracted in LALBarycenter.c
    }
    else{
      GSL_SET_COMPLEX(&emitdt, emit.deltaT, 0.);
    }

    // fill in the output training buffer
    gsl_vector_complex_set(wv, i, emitdt);
  }
}


void Binary_Barycenter_Waveform(gsl_vector_complex *wv,
                                const gsl_vector *timestamps,
                                const double *params,
                                const std::string model_tag){
  int n = timestamps->size;
  REAL8 w0 = params[0];  // angle of periastron (0 -> 2*pi)
  REAL8 T0 = params[1];  // time of periastron  (0 -> 1) - re-scale to Pb range
  REAL8 ecc = params[2]; // eccentricity

  // period is just set to be the whole timespan
  REAL8 Pb = gsl_vector_get(timestamps, timestamps->size-1) - gsl_vector_get(timestamps, 0);

  T0 *= Pb; // rescale to be between 0 and Pb

  // set asini to be 1
  REAL8 asini = 1.;

  BinaryPulsarInput binInput;
  BinaryPulsarOutput binOutput;

  char model[256] = "BT";
  
  // variables for calculating barycenter time delay
  PulsarParameters *pars = (PulsarParameters*)XLALMalloc(sizeof(PulsarParameters *));

  PulsarAddParam( pars, "OM", &w0, PULSARTYPE_REAL8Vector_t );
  PulsarAddParam( pars, "T0", &T0, PULSARTYPE_REAL8Vector_t );
  PulsarAddParam( pars, "ECC", &ecc, PULSARTYPE_REAL8Vector_t );
  PulsarAddParam( pars, "A1", &asini, PULSARTYPE_REAL8Vector_t );
  PulsarAddParam( pars, "PB", &Pb, PULSARTYPE_REAL8Vector_t );
  PulsarAddParam( pars, "BINARY", model, PULSARTYPE_string_t );

  for ( int i=0; i < n; i++ ){
    binInput.tb = gsl_vector_get(timestamps, i);

    // calculate binary time delay
    XLALBinaryPulsarDeltaTNew( &binOutput, &binInput, pars );

    // fill in the output training buffer
    gsl_complex emitdt;
    GSL_SET_COMPLEX(&emitdt, binOutput.deltaT, 0.);
    gsl_vector_complex_set(wv, i, emitdt);
  }
}
