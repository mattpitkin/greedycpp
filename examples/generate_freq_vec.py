#! /usr/bin/env python

import sys
import numpy
from optparse import OptionParser

# kludge to simulate importing from lal
class lal:
    LAL_C_SI = 299792458
    LAL_G_SI = 6.67259e-11
    LAL_MSUN_SI = 1.98892e+30
    LAL_MTSUN_SI = 4.925492321898864e-06

#
#   Functions to estimate PN expansion
#
def schwarzschild_fisco(m1, m2):
    return lal.LAL_C_SI**3. / (6.**(3./2) * numpy.pi * (m1+m2)*lal.LAL_MSUN_SI
            * lal.LAL_G_SI)

# The following are from Poisson & Will
def eta_from_m1m2(m1, m2):
    return m1*m2 / (m1+m2)**2.

def mchirp_from_m1m2(m1, m2):
    return eta_from_m1m2(m1, m2)**(3./5)*(m1+m2)

def so_coupling(m1, m2, s1z, s2z):
    m1 = m1 * lal.LAL_MTSUN_SI
    m2 = m2 * lal.LAL_MTSUN_SI
    eta = eta_from_m1m2(m1, m2)
    M = m1+m2
    return sum([(113 * (m_i/M)**2. + 75 * eta)*s_i for m_i, s_i in
        [(m1, s1z), (m2, s2z)]])/12.

def ss_coupling(m1, m2, s1, s2):
    eta = eta_from_m1m2(m1*lal.LAL_MTSUN_SI, m2*lal.LAL_MTSUN_SI)
    return (eta/48.) * (-247 * numpy.dot(s1, s2) + 721 * s1[2] * s2[2])

def t0PN(m1, m2, f):
    """
    Gives the time-to-coalesence using Newtonian estimate. Also known as tau0.

    Parameters
    ----------
    m1: float
        Mass of the larger object in solar masses
    m2: float
        Mass of the smaller object in solar masses
    f: float
        Starting frequency in Hz

    Returns
    -------
    tau0: float
        The time-to-coalesence from the starting frequency, in seconds.
    """
    mchirp = mchirp_from_m1m2(m1, m2) * lal.LAL_MTSUN_SI
    return (5./256)* mchirp * (numpy.pi * mchirp * f)**(-8./3)

def t1PN(m1, m2, f):
    """
    Gives the 1PN correction of the time-to-coalesence.

    Parameters
    ----------
    m1: float
        Mass of the larger object in solar masses
    m2: float
        Mass of the smaller object in solar masses
    f: float
        Starting frequency in Hz

    Returns
    -------
    t1PN: float
        The 1PN correction of the time-to-coalesence in seconds
    """
    mchirp = mchirp_from_m1m2(m1, m2) * lal.LAL_MTSUN_SI
    eta = eta_from_m1m2(m1, m2)
    M = (m1+m2) * lal.LAL_MTSUN_SI
    return  (4./3) * (743./336 + (11./4)*eta) * (numpy.pi * M * f)**(2./3)

def t1_5PN(m1, m2, s1, s2, f):
    """
    Gives the 1.5PN correction of the time-to-coalesence.

    Parameters
    ----------
    m1: float
        Mass of the larger object in solar masses
    m2: float
        Mass of the smaller object in solar masses
    s1: array of floats
        Spin components of the larger object
    s2: array of floats
        Spin components of the smaller object
    f: float
        Starting frequency in Hz

    Returns
    -------
    t1_5PN: float
        The 1.5PN correction of the time-to-coalesence in seconds
    """
    M = (m1+m2) * lal.LAL_MTSUN_SI
    beta = so_coupling(m1, m2, s1[2], s2[2])
    return -(8./5) * (4*numpy.pi - beta) * (numpy.pi * M * f)

def t2PN(m1, m2, s1, s2, f):
    """
    Gives the 2PN correction of the time-to-coalesence.

    Parameters
    ----------
    m1: float
        Mass of the larger object in solar masses
    m2: float
        Mass of the smaller object in solar masses
    s1: array of floats
        Spin components of the larger object
    s2: array of floats
        Spin components of the smaller object
    f: float
        Starting frequency in Hz

    Returns
    -------
    t2PN: float
        The 2PN correction of the time-to-coalesence in seconds
    """
    eta = eta_from_m1m2(m1, m2)
    sigma = ss_coupling(m1, m2, s1, s2)
    M = (m1+m2) * lal.LAL_MTSUN_SI
    return 2*(3058673./1016064 + (5429./1008)*eta +
        (617./144)*eta**2. - sigma) * (numpy.pi * M * f)**(4./3)

def t_of_F(m1, m2, s1, s2, f, tc=0., order=4):
    """
    Gives the time-to-coalesence at the given frequency using the PN
    approximation to the given order.

    Parameters
    ----------
    m1: float
        Mass of the larger object in solar masses
    m2: float
        Mass of the smaller object in solar masses
    s1: array of floats
        Spin components of the larger object
    s2: array of floats
        Spin components of the smaller object
    f: float
        Starting frequency in Hz
    tc: float
        Time at coalesence in seconds. Default is 0.
    order: int
        Order of largest correction to use (xPN = order/2). Currently supports
        up to 4 (= 2PN).

    Returns
    -------
    tc - t: float
        The amount of time from the given frequency until coalesence,
        in seconds.
    """
    PNs = numpy.array([
            1., 0., t1PN(m1, m2, f), t1_5PN(m1, m2, s1, s2, f),
        t2PN(m1, m2, s1, s2, f)])
    if order > (len(PNs)-1):
        raise ValueError("order must be <= %i" %(len(PNs)-1))
    return tc - t0PN(m1, m2, f)*PNs[:order+1].sum()

def estimate_duration(m1, m2, s1, s2, f0, f, order=4):
    """
    Estimates the inspiral duration between an initial frequency f0 and a final
    frequency f using post-Newtonian approximation to the given order.

    Parameters
    ----------
    m1: float
        Mass of the larger object in solar masses
    m2: float
        Mass of the smaller object in solar masses
    s1: array of floats
        Spin components of the larger object
    s2: array of floats
        Spin components of the smaller object
    f0: float
        Starting frequency in Hz
    f: float
        Ending frequency in Hz
    order: int
        Order of largest correction to use (xPN = order/2). Currently supports
        up to 4 (= 2PN).

    Returns
    -------
    duration: float
        An estimate of the duration of an inspiral waveform between f0 and f,
        in seconds.
    """
    # we do t(f0) - t(f) since t_of_F with tc = 0 returns the number
    # of seconds before coalesence
    return t_of_F(m1, m2, s1, s2, f, order=order) - \
        t_of_F(m1, m2, s1, s2, f0, order=order)


def get_seglength(dur):
    """
    Given a non-integer duration, finds the segment length that is the closest
    power of 2 that will contain the duration.
    """
    return 2**numpy.ceil(numpy.log2(dur))


#
#   Main
#
parser = OptionParser()
parser.add_option('-b', '--bank-points', help='Required. Text file containing mass points')
parser.add_option('-o', '--output-file', help='Required. Name of text file to write the frequency vector out to.')
parser.add_option('-f', '--f-low', type=float, help='Required. Starting frequency to use')
parser.add_option('-F', '--f-stop', type=float, default=None, help='Optional. Ending frequency to use. If None specified, the largest Schwarzschild ISCO in the bank will be used.')
parser.add_option('-S', '--default-segment-length', type=float, default=None, help='Optional. The inverse of this gives the default frequency step to use. That is, the frequency vector will be constructed such that every point is a multiple of the default frequency step. If not specified, the longest duration in the bank rounded up to the nearest power of two will be used.')

opts, _ = parser.parse_args()

# check opts
if opts.bank_points is None:
    raise ValueError("Must specify bank-points file")
if opts.output_file is None:
    raise ValueError("Must specify output-file")
if opts.f_low is None:
    raise ValueError("Must specify f-low")

# load the bank points
print >> sys.stdout, "Loading bank..."
bank = numpy.loadtxt(opts.bank_points)
# FIXME: eventually add spins to bank files
s1 = [0., 0., 0.]
s2 = [0., 0., 0.]

# if no ending specified, finding highest isco
if opts.f_stop is None:
    opts.f_stop = schwarzschild_fisco(bank[:,0], bank[:,1]).max()
    print >> sys.stdout, "Using %f as stopping frequency..." % opts.f_stop

####
# FIXME: Assuming that can just use the longest waveform
# should really check this
durs = numpy.array([estimate_duration(m1, m2, s1, s2, opts.f_low,
            opts.f_stop) for m1, m2 in bank])
use_idx = durs.argmax()
m1, m2 = bank[use_idx,:]
####

if opts.default_segment_length is not None:
    default_seg_len = opts.default_segment_length
else:
    default_seg_len = get_seglength(durs[use_idx])
    print "Using %is as the default segment length." % default_seg_len
default_df = 1./default_seg_len

# now sweep through the frequency band, adjusting the step as needed
print >> sys.stdout, "Calculating frequency vector..."
freq_vec = []
freq = opts.f_low
while freq < opts.f_stop:
    freq_vec.append(freq)
    # figure out the next df to use:
    # this is 1 / the segment length needed to contain the waveform with the
    # longest duration from the current frequency to the terminating frequency
    # in the bank
    #next_df = 1./get_seglength(max([estimate_duration(m1, m2, s1, s2, freq,
    #        opts.f_stop) for m1, m2 in bank]))
    next_df = numpy.floor(default_seg_len/estimate_duration(m1, m2, s1, s2, freq,
        opts.f_stop)) * default_df
    freq += next_df
    print >> sys.stdout, "%f %e\r" %(freq, next_df),
    sys.stdout.flush()
print ""
# add the last point for full coverage
freq_vec.append(opts.f_stop)

# write the frequency vector out and exit
print >> sys.stdout, "Saving..."
numpy.savetxt(opts.output_file, numpy.array(freq_vec))

N = len(freq_vec)
Norig = (opts.f_stop - opts.f_low) / default_df
print >> sys.stdout, "Number of frequency points: %i" % N
print >> sys.stdout, "Number without dynamic adjustment: %i" % Norig
print >> sys.stdout, "Savings factor: %f" %(Norig/N)
