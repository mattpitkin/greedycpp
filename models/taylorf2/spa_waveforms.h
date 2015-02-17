//#include <iostream>
#include <stdio.h>
#include <math.h>

/*----------------------------------------------------------------------------
 3.5PN inspiral waveform (TaylorF2) in the stationary phase approximation
 Used in both ROQ and Full computations
 -----------------------------------------------------------------------------*/

void TF2_Waveform(double &TF2_Real, double &TF2_Imag, double *theta, double frq, double amp, double PN)
{
    const double gammaE = 0.5772156649015328; // Euler–Mascheroni constant

    double Mtot    = theta[0] + theta[1];
    double Msym    = theta[0]*theta[1]/(Mtot*Mtot);
    double Mchirp  = pow(theta[0]*theta[1], 3./5.)*pow(Mtot, -.1/5.);
    
    double Vchar;                      // Characteristic velocity
    double v2,v3,v4,v5,v6,v7;          // Powers of Characteristic velocity
    double TF2_Phase;                  // PN GW phase
    double Msym2, Msym3;               // Powers of the symmetric mass
    double cA, q2, q4,q5,q6,q7;
    
    //printf("%4.6e %4.6e  %4.6e  %4.6e %4.6e\n\n", theta[0],theta[1], theta[2],theta[3],amp);//exit(0);
    
    /*---- Powers of the Frequency and symmetric mass---*/
    Vchar = pow(M_PI*Mtot*frq, 1./3.);
    v2 = Vchar*Vchar;
    v3 = v2*Vchar;
    v4 = v2*v2;
    v5 = v4*Vchar;
    v6 = v3*v3;
    v7 = v6*Vchar;
    
    Msym2 = Msym*Msym;
    Msym3 = Msym*Msym2;
    
    /*---- PN Phase coefficients ----*/
    cA = 3./(128.*Msym*v5);

    q2 = 20.*(743./336. + 11./4.*Msym)/9.;
    q4 = 10*(3058673./1016064. + 5429./1008.*Msym + 617./144.*Msym2);
    q5 = M_PI*( 38645./756. - 65./9.*Msym)*(1. + 3.*log(Vchar));
    
    q6 = 11583231236531./4694215680. - 640./3.*(M_PI*M_PI)
    - 6848./21.*(gammaE + log(4.*Vchar))+ (2255.*(M_PI*M_PI)/12.-15737765635.0/3048192.0 )*Msym
    + 76055./1728.*Msym2 - 127825./1296.*Msym3;
    
    q7 = M_PI*(77096675./254016. + 378515./1512.*Msym - 74045./756.*Msym2);
    
    if(PN==0)
        TF2_Phase = 2.*M_PI*frq*theta[2] - theta[3] - M_PI/4. + cA; // 0 order PN chirp waveform
    else if(PN==3.5)
        TF2_Phase = 2.*M_PI*frq*theta[2] - theta[3] - M_PI/4. + cA*(1. + q2*v2 - 16.*M_PI*v3 + q4*v4 + q5*v5 +q6*v6 +q7*v7);
    else
    {
        std::cerr << "Choice not coded yet!" << std::endl;
        exit(0);
    }
    
    TF2_Real = amp*pow(frq, -7./6.)*cos(TF2_Phase);
    TF2_Imag = amp*pow(frq, -7./6.)*sin(TF2_Phase);
    
   // printf(" %4.6e %4.6e\n", TF2_Real,TF2_Imag);//exit(0);
    
}

// this routine is interface with greedy routine -- returns gsl data type //
void TF2_FullWaveform(gsl_vector_complex *wv, double *params, const gsl_vector *xQuad, double amp, double PN)
{

    double *more_params;
    more_params = new double[4]; // (m1,m2,tc,phi_c)

    more_params[0] = params[0];
    more_params[1] = params[1];
    more_params[2] = 0.0;  // dummy variable (tc in waveform generation)
    more_params[3] = 0.0;  // dummy variable (phi_c in waveform generation)


    // parameter list such that (m1(param),m2(param)) is a unique point in parameter space
    double TS_r = 0.0;
    double TS_i = 0.0;
    gsl_complex zM;

    for(int cols = 0; cols < xQuad->size; cols++)
    {
        TF2_Waveform(TS_r, TS_i, more_params, xQuad->data[cols], amp, PN);
        GSL_SET_COMPLEX(&zM, TS_r, TS_i);
        gsl_vector_complex_set(wv,cols,zM);
    }

    delete[] more_params;
}


