"""
Module to do sanity tests of output directory...

1) basis orthogonality test
2) basis accuracy test (training space file needs to be present)


Can run at command line from an application directory, e.g.:
  $ python sanity_tests.py OUTDIR1
"""

# TODO: write B-matrix style and roq sanity test

import os, re, sys
import numpy as np
import greedy_utils

def ortho_test(outdir,xmin=None,xmax=None,quad_points=None):

  B, weights, nodes = greedy_utils.load_info(outdir,xmin,xmax,quad_points)

  result = np.dot(weights * B.T.conj(),B) 
  err = result - np.eye(result.shape[1])
  print "othogonality error %1.15e" % np.max(np.abs(err))

def basis_accuracy_test(outdir):

  B, weights, nodes = greedy_utils.load_info(outdir)
  TS = greedy_utils.load_ts(outdir)
  quad_points, evaluations = TS.real.shape 

  err_app = np.zeros(evaluations)

  # TODO: shouldn't use for-loop here
  for ii in range(evaluations):
    proj_coeff = np.dot(weights * B.T.conj(),TS[:,ii])
    err_h = TS[:,ii] - np.dot(B,proj_coeff) # not conj is good! we just want the sum (no weights either)
    err_app[ii] = np.sqrt( np.dot( weights * err_h,err_h.conj()) ).real

  import matplotlib.pyplot as plt
  plt.semilogy(range(evaluations),err_app,'ro',label='rb error')
  #plt.show()

  ## to show other info ##
  #import matplotlib.pyplot as plt
  #plt.plot(nodes,TS[:,1])
  #plt.show()
  #plt.plot(nodes,err_h)
  #plt.show()

def eim_accuracy_test(outdir):

  B, weights, nodes = greedy_utils.load_info(outdir)
  TS = greedy_utils.load_ts(outdir)
  quad_points, evaluations = TS.real.shape
  invV, eim_indx, eim_nodes = greedy_utils.load_eim(outdir)

  err_app = np.zeros(evaluations)

  ## build the empirical interpolant ## 
  coeffs = np.dot(invV, TS[eim_indx,:])
  ts_eim = np.dot(B,coeffs)

  # TODO: shouldn't use for-loop here
  for ii in range(evaluations):
    err_h = TS[:,ii] - ts_eim[:,ii]
    err_app[ii] = np.sqrt( np.dot(weights * err_h, err_h.conj() ) ).real

  import matplotlib.pyplot as plt
  plt.semilogy(range(evaluations),err_app,'bo',label='eim error')
  plt.xlabel('ts index')
  plt.ylabel('error')
  plt.legend()
  plt.show()

if __name__=="__main__":

  try:
    outdir = sys.argv[1]
  except:
    raise Exception("Must specify an output directory")

  ortho_test(outdir)
  basis_accuracy_test(outdir)
  eim_accuracy_test(outdir)

