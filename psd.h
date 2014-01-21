//#include <iostream>
#include <stdio.h>
#include <math.h>

/*--- POWER SPECTRAL DENSITY ---*/
void ComputePSD(double *Sn, int noise_flag, double f_low, double deltaF, int FDlen)
{

    double f,x,fo,x2,x4, poly; // Reference frequency in Hz

    for (int i=0;i<FDlen;i++)
    {
        f = i*deltaF + f_low;
        switch (noise_flag)
        {
            case 1: // White noise
                Sn[i]=1.;
                break;
            case 2: // LIGO noise curve
                fo = 150;
                x = f/fo;
                x2 = x*x;
                Sn[i] = 9E-46*(pow(4.49*x,-56) + 0.16*pow(x,-4.52)+ 0.52 + 0.32*x2);
                break;
            case 3: //AdLIGO noise curve
                fo = 215;
                x = f/fo;
                x2 = x*x;
                x4 = x2*x2;
                poly = 1. - x2 + x4*0.5;
                poly /= 1 + x2*0.5;
                Sn[i] = 1E-49*(pow(x,-4.14) - 5/x2+ 111*poly);
                break;
            default:
                std::cerr << "Noise choice not coded" << std::endl;
                exit(0);
        }

        //fprintf(stdout,"noise at f = %f is Sn = %e\n",f,Sn[i]);
    }

}

