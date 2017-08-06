#!/usr/bin/env python

"""
Generate training points for barycentering: right ascension and declination pairs
drawn randomly across the sky
"""

import numpy as np
import argparse

# main function
if __name__=='__main__':
  description = """This script will generate a set of training points in w0, T0, and eccentricity"""

  parser = argparse.ArgumentParser( description = description )
  parser.add_argument("-N", "--ntraining", dest="ntraining", default=10000, type=int, help="Number of training samples")
  parser.add_argument("-o", "--outfile", dest="outfile", default="TS_points.txt", help="Output training set file name")

  # parse input options
  opts = parser.parse_args()

  # generate training points
  w0 = 2.*np.pi*np.random.rand(opts.ntraining) # uniform from 0 -> 2pi
  T0 = np.random.rand(opts.ntraining)          # uniform from 0 -> 1
  ecc = 0.5*np.random.rand(opts.ntraining)     # uniform from 0 -> 0.5

  fp = open(opts.outfile, 'w')
  for i in range(opts.ntraining):
    fp.write('%.14f\t%.14f\t%.14f\n' % (w0[i], T0[i], ecc[i]))
  fp.close()
