#! /usr/bin/env python

import numpy
import lalsimulation as lalsim

psd_model = 'iLIGOSRD'
use_freq_vec = True
if use_freq_vec:
    freq_points = numpy.loadtxt('frequency_vector.txt')
else:
    flow = 40.
    fhigh = 1024.
    df = 1./64
    N = int((fhigh - flow)/df)
    freq_points = numpy.arange(N)*df + flow

f = getattr(lalsim, 'SimNoisePSD%s' % psd_model)
psd = map(f, freq_points)
    
print "Num freq points: %f" % len(freq_points)

numpy.savetxt('ASD.txt', numpy.sqrt(psd))
