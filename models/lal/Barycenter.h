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
  REAL8 T0 = params[0];    // time of periastron  (0 -> 1) - re-scale to Pb range
  REAL8 ecc = params[1];   // eccentricity

  // period is just set to be the whole timespan
  REAL8 Pb = gsl_vector_get(timestamps, timestamps->size-1) - gsl_vector_get(timestamps, 0);

  T0 *= Pb; // rescale to be between 0 and Pb

  // deduce whether wanting to output the time delay derivative
  int sinorcos = lal_help::get_binary_barycenter_tags(model_tag);

  REAL8 phase = 0., orbits = 0., U = 0.;
  INT4 norbits = 0;

  for ( int i=0; i < n; i++ ){
    orbits = (gsl_vector_get(timestamps, i) - T0)/Pb;
    norbits = (INT4)orbits;
    if ( orbits < 0. ) norbits--;
    phase = LAL_TWOPI*(orbits - (REAL8)norbits);
    XLALComputeEccentricAnomaly( phase, ecc, &U );

    // fill in the output training buffer
    gsl_complex emitdt;

    if ( sinorcos == 1 ){ GSL_SET_COMPLEX(&emitdt, sin(U), 0.); }
    else { GSL_SET_COMPLEX(&emitdt, cos(U), 0.); }

    gsl_vector_complex_set(wv, i, emitdt);
  }
}

