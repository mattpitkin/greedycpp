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
#include <time.h>

#include <lal/LALConstants.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/SFTutils.h>
#include <lal/BinaryPulsarTiming.h>
#include <lal/Date.h>

#include <omp.h>

#ifdef MODEL_TEMPO
#include <tempo2.h>
#endif

// global variables for barycentering (so ephemeris files are only read in once)
EphemerisData *edat = NULL;
TimeCorrectionData *tdat = NULL;
int usetempo = 0;
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
  
#ifdef MODEL_TEMPO
  struct pulsar *psr = NULL;
#endif

  if ( !edat && !tdat && !usetempo ){
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

#ifdef MODEL_TEMPO
    // check whether using LAL or tempo2...
    if ( !strcmp("TEMPO", vals[4].c_str()) ){
      usetempo = 1;
      fprintf(stderr, "using TEMPO model\n");
    }
#else
    usetempo = 0;
#endif

    if ( !usetempo ){
      // set earth and sun ephemeris files
      char earthfile[256], sunfile[256];
      snprintf(earthfile, sizeof(char)*256, "earth00-19-%s.dat.gz", vals[1].c_str());
      snprintf(sunfile, sizeof(char)*256, "sun00-19-%s.dat.gz", vals[1].c_str());
      edat = XLALInitBarycenter( earthfile, sunfile );
    }

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
    if ( !usetempo ) { tdat = XLALInitTimeCorrections( timecorrfile ); }

    if ( !strcmp("NOSHAPIRO", vals[3].c_str()) ){
      noshapiro = 1; // do not include Shapiro delay in barycenter calculation
    }
    #ifdef USE_OPENMP
    }
    #pragma omp barrier
    #endif
  }

#ifdef MODEL_TEMPO
  // setup pulsar
  if ( usetempo ){
    std::vector<std::string> vals = lal_help::get_barycenter_tags(model_tag);

    // initialise pulsar
    psr = (pulsar *)malloc(sizeof(pulsar)*1);
    MAX_OBSN = n; // define MAX_OBSN (in tempo2.h) as the number of observations
    initialiseOne(psr, 1, 1); // initialise pulsar
  
    strcpy(psr[0].JPL_EPHEMERIS, getenv(TEMPO2_ENVIRON));
    char epfile[256];
    snprintf(epfile, sizeof(char)*256, "/ephemeris/%s.1950.2050", vals[1].c_str());
    strcpy(psr[0].ephemeris, vals[1].c_str());
    strcat(psr[0].JPL_EPHEMERIS, epfile);

    if ( !strcmp("TCB", vals[2].c_str()) ){ psr[0].units = SI_UNITS; }
    if ( !strcmp("TDB", vals[2].c_str()) ){ psr[0].units = TDB_UNITS; }

    // set the site (assume that LIGO sites have been added to the TEMPO2 observatories file)
    if ( strstr(vals[0].c_str(), "H1") != NULL ){
      strcpy(psr[0].obsn[0].telID, "HANFORD");
    }
    else if ( strstr(vals[0].c_str(), "L1") != NULL ){
      strcpy(psr[0].obsn[0].telID, "LIVINGSTON");
    }
    else if ( strstr(vals[0].c_str(), "V1") != NULL ){ 
      strcpy(psr[0].obsn[0].telID, "VIRGO");
    }
    else{
      fprintf(stderr, "Error... detector not found.\n");
      exit(-1);
    }
  }
#endif

  if ( !usetempo ){
    // variables for calculating barycenter time delay
    EarthState earth;
    EmissionTime emit;
    BarycenterInput baryinput;

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

      //if ( i == 523 ){ fprintf(stderr, "time = %.16lf\n", GSL_REAL(emitdt)); }
      // fill in the output training buffer
      gsl_vector_complex_set(wv, i, emitdt);
    }
  }
#ifdef MODEL_TEMPO
  else{ // use tempo2 for barycentring
    REAL8 batdt = 0.;
    
    // set pulsar position
    psr[0].param[param_raj].val[0] = ra;
    psr[0].param[param_decj].val[0] = dec;
    psr[0].param[param_px].val[0] = 0.; // no parallax

    psr[0].correctTroposphere = 0;
    vectorPulsar(psr, 1);

    psr[0].nobs = n;
    for ( int i=0; i < n; i++ ){
      psr[0].obsn[i].delayCorr = 1;
      psr[0].obsn[i].clockCorr = 1;

      // set SAT value of observation
      REAL8 thistime = gsl_vector_get(timestamps, i);
      REAL8 fracsec = thistime - floor(thistime);     // if not an integer get remaining fraction of seconds
      REAL8 mjd = 0.; // time as Modified Julian Date
      
      // this time is a GPS time, so this needs to be converted to UTC and then to an MJD
      // NOTE: this requires the times to be integer seconds
      struct tm utc;
      XLALGPSToUTC( &utc, (INT4)floor(thistime) ); // convert GPS to UTC
      mjd = XLALConvertCivilTimeToMJD( &utc );     // convert UTC into MJD format
      mjd += fracsec/86400.;                       // add on fractional seconds
  
      psr[0].obsn[i].sat = (long double)mjd;

      if ( i > 0 ){
        strcpy(psr[0].obsn[i].telID, psr[0].obsn[0].telID);
      }
    }

    formBatsAll(psr, 1); // should contain everything that's required and in the right order

    for ( int i=0; i < n; i++ ){
      gsl_complex emitdt;
      REAL8 shapirodelay = 0.;

      // get barycenter time delat
      batdt = psr[0].obsn[i].correctionTT_TB + psr[0].obsn[i].roemer;
      
      if ( !noshapiro ){ // include Shapiro delay
        shapirodelay = psr[0].obsn[i].shapiroDelaySun; // only include Sun
      }

      //fprintf(stderr, "correction = %.16lf  ", getCorrectionTT(psr[0].obsn));
      GSL_SET_COMPLEX(&emitdt, batdt - shapirodelay, 0.); // subtract shapiro as in TEMPO

      // fill in the output training buffer
      //if ( i == 523 ){ fprintf(stderr, "time = %.16lf\n", GSL_REAL(emitdt)); }
      gsl_vector_complex_set(wv, i, emitdt);
    }
    destroyOne(psr); // free memory
  }
#endif
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

