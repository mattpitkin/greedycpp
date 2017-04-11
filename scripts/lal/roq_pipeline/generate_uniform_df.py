#!/usr/bin/env python

import lal
import lalsimulation as LS
import numpy as np


def Mfun(Mc, eta):
    return Mc*eta**(-3.0/5.0)


def qfun(eta):
    return (1.0 + np.sqrt(1.0 - 4.0*eta) - 2.0*eta) / (2.0*eta)


def m1fun(M, q):
    return M*1.0/(1.0+q)


def m2fun(M, q):
    return M*q/(1.0+q)


def ComputeDeltaF(Mc, eta, chi=0.99, f_min=20):
    Mtot = Mfun(Mc, eta)
    q = qfun(eta)
    # print Mtot, q
    m1 = m1fun(Mtot, q) * lal.MSUN_SI
    m2 = m2fun(Mtot, q) * lal.MSUN_SI

    T = LS.SimIMRSEOBNRv2ROMDoubleSpinHITimeOfFrequency(f_min, m1, m2,
                                                        chi, chi)
    # print 'T = 1/df[Hz]', T
    # print 'df [Hz]', 1.0 / T
    return round(T)

Mc_mins = [1.4, 2.2, 3.4, 5.2, 7.9, 12.3]
dfs = [ComputeDeltaF(Mc, 0.09, chi=0.99, f_min=20) for Mc in Mc_mins]
print Mc_mins
print dfs
