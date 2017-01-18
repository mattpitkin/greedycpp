#!/usr/bin/env python

from generate_trainingpoints import generate_sampling, parse_param_list
import generate_param_samples as gps

import numpy as np
import random, argparse
#from optparse import OptionParser

# Parameter transformation functions
def etafun(q):
  return q/(1.0 + q)**2
def qfun(eta):
  return (1.0 + np.sqrt(1.0 - 4.0*eta) - 2.0*eta) / (2.0*eta)
def m1fun(M,q):
  return M*1.0/(1.0+q)
def m2fun(M,q):
  return M*q/(1.0+q)
def Mchirpfun(M, eta):
  return M*eta**(3.0/5.0)
def Mfun(Mc, eta):
  return Mc*eta**(-3.0/5.0)
def Mcnu_to_m1m2(Mc, eta):
  M = Mfun(Mc, eta)
  q = qfun(eta)
  return M*q/(1.0+q), M/(1.0+q)

def main():

  parser=argparse.ArgumentParser(description='seob sampler',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument("-i","--input",dest="input", action="store", type=str,help="input file containing sampling ranges, modes and types")
  parser.add_argument("-o","--out",dest="fout", action="store", default="./trainingpoints.txt", type=str,help="output file to store the trainingpoints")
  parser.add_argument('--sampler',dest='sampler',type=str,default='file',choices=['random','file'],
                      help='strategy for sampling the parameter space. If file, use the input file. If random, N is the TOTAL number of samples and ignore full vs no-full.')
  parser.add_argument('--full', dest='full',
                      help='TS over full parameter domain (full=True). Ignored if sampler = random.',action='store_true')
  parser.add_argument('--no-full', dest='full',
                     help='TS over boundary of parameter domain (full=False). Ignored if sampler = random.',action='store_false')
  parser.add_argument('--N', dest='N',type=int,default=40,
                      help='number of RANDOM points per parameter space dimension')

  #(opts,args) = parser.parse_args()
  opts = parser.parse_args()

  full  = opts.full
  f_in  = opts.input
  f_out = opts.fout
  N     = opts.N
  sampler = opts.sampler

  print full
  print f_in
  print f_out
  print N
  print sampler


  if sampler == 'file':
    ts_size, ts = generate_sampling('junk.txt',param_list=f_in,file_to_extend="",vector_files=[],return_ts=True)
    print "Trainingpoints written to \"%s\""%f_out
    rows = ts.shape[0]
    cols = ts.shape[1]

    ts_mins = [ts.min(axis=0)[i] for i in np.arange(cols)]
    ts_maxs = [ts.max(axis=0)[i] for i in np.arange(cols)]

    # sanity check and direct bounds checking
    intervals = np.genfromtxt(f_in,names=True)
    for i, name in enumerate(intervals.dtype.names):
      print("Parameter %s, min = %f, max = %f"%(name,intervals[name][0],intervals[name][1]))
      print("Parameter %s, min(TS) = %f, max(TS) = %f"%(name,ts_mins[i],ts_maxs[i]))

    def on_boundary(x):
      """ is the point x on the training set's boundary?"""
      bdry = False
      for i, x_i in enumerate(x):
        xi_min = ts_mins[i]
        xi_max = ts_maxs[i]
        #print("i = %i with xi =%f and max/min = %f, %f"%(i,x_i,xi_max,xi_min))
        if np.abs(x_i-xi_min)<1.e-5 or np.abs(x_i-xi_max)<1.e-5 :
          bdry = True
          break
      return bdry

    # (eta, Mc, chi1, chi2) -> (m1, m2, chi1, chi2) with additional restrictions
    Mmax = 0.0
    f = open(f_out, 'w')

    for i in np.arange(rows):
      x = ts[i,:]
      Mc  = x[0]
      eta = x[1]
      chi1 = x[2]
      chi2 = x[3]
      m1, m2 = Mcnu_to_m1m2(Mc, eta)

      if m1 >= 1. and m2 >= 1.:
        mm = m1+m2
        if mm > Mmax:
          Mmax = mm  # update maximum total mass
        if full:
          f.write('%f %f %f %f\n'%(m1, m2, chi1, chi2))
        else: #Only write point if we are on a boundary surface of the domain
          if on_boundary(x):
            f.write('%f %f %f %f\n'%(m1, m2, chi1, chi2))
    f.close()
    print 'Mmax', Mmax

  elif sampler == 'random':

    param_info, params_low, params_high, params_num, params_type = parse_param_list(f_in)
    a_ranges = np.array([params_low,params_high],dtype=float)

    print params_low
    print params_high
    print a_ranges

    if f_out.rfind('/') == -1:
      postfix = f_out
    else:
      postfix = f_out[f_out.rfind('/')+1:] # postfix on file needed if multiple files are being created at the same time through other validations
    tmp_filename = 'tmp_McEtaSampling.'+postfix
    print "tmp_filename = %s"%tmp_filename
    random_pts = gps.generate_sampling(tmp_filename,a_ranges,total_picks=N,return_picks=True)

    rows      = random_pts.shape[0]
    param_dim = random_pts.shape[1]
    print rows
    print param_dim

    # (eta, Mc, chi1, chi2) -> (m1, m2, chi1, chi2) with additional restrictions
    Mmax = 0.0
    f = open(f_out, 'w')

    for i in np.arange(rows):
      x   = random_pts[i,:]
      Mc  = x[0]
      eta = x[1]
      chi1 = x[2]
      chi2 = x[3]
      m1, m2 = Mcnu_to_m1m2(Mc, eta)

      if m1 >= 1. and m2 >= 1.:
        mm = m1+m2
        if mm > Mmax:
          Mmax = mm  # update maximum total mass
        f.write('%f %f %f %f\n'%(m1, m2, chi1, chi2))
        #if full:
        #  f.write('%f %f %f %f\n'%(m1, m2, chi1, chi2))
        #else: #Only write point if we are on a boundary surface of the domain
        #  if on_boundary(x):
        #    f.write('%f %f %f %f\n'%(m1, m2, chi1, chi2))
    f.close()
    print 'Mmax', Mmax

if __name__ == "__main__":
	# START THE MAIN FUNCTION IF RUN AS A SCRIPT. OTHERWISE, JUST LOADS THE CLASS
	# AND FUNCTION DEFINITIONS
	exit(main())
