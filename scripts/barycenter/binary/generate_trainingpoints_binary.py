#!/usr/bin/env python

"""
Generate training points for binary barycentering: time of periastron and eccentricity pairs
drawn randomly
"""

import numpy as np
import argparse

# main function
if __name__=='__main__':
    description = """This script will generate a set of training points in T0 and eccentricity"""

    parser = argparse.ArgumentParser( description = description )
    parser.add_argument("-N", "--ntraining", dest="ntraining", default=10000, type=int, help="Number of training samples")
    parser.add_argument("-e", "--max-ecc", dest="maxecc", default=0.1, type=float, help="Maximum eccentricity of training samples")
    parser.add_argument("-o", "--outfile", dest="outfile", default="TS_points.txt", help="Output training set file name")

    # parse input options
    opts = parser.parse_args()

    # generate training points
    T0 = np.random.rand(opts.ntraining)              # uniform from 0 -> 1
    ecc = opts.maxecc*np.random.rand(opts.ntraining) # uniform from 0 -> max eccentricity

    fp = open(opts.outfile, 'w')
    for i in range(opts.ntraining):
        fp.write('%.14f\t%.14f\n' % (T0[i], ecc[i]))
    fp.close()
