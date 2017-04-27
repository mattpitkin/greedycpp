#!/usr/bin/env python

"""
Generate training points for barycentering: right ascension and declination pairs
drawn randomly across the sky
"""

import numpy as np
import argparse

# main function
if __name__=='__main__':
  description = """This script will generate a set of training points in right ascension and declination"""

  parser = argparse.ArgumentParser( description = description )
  parser.add_argument("-N", "--ntraining", dest="ntraining", default=10000, type=int, help="Number of training samples")
  parser.add_argument("-o", "--outfile", dest="outfile", default="TS_points.txt", help="Output training set file name")

  # parse input options
  opts = parser.parse_args()

  # generate training points (uniformly from across the sky)
  ra = 2.*np.pi*np.random.rand(opts.ntraining)
  dec = -(np.pi/2.) + np.arccos(2.*np.random.rand(opts.ntraining) - 1.)

  fp = open(opts.outfile, 'w')
  for i in range(opts.ntraining):
    fp.write('%.14f\t%.14f\n' % (ra[i], dec[i]))
  fp.close()
