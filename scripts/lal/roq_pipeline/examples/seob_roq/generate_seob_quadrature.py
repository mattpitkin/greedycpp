#!/usr/bin/env python

import numpy as np
import argparse
from generate_quadrature_rule import GravitationalWaveAdaptive
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

  parser  = argparse.ArgumentParser(description='GW quadrature rules',\
           formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument("-i","--input",dest="input", action="store", nargs='+', type=str, help="input files containing sampling ranges")
  parser.add_argument("-l","--flow",dest="flow", action="store", type=float,help="lower frequency")
  parser.add_argument("-u","--fup",dest="fup", action="store", type=float,help="lower frequency")
  parser.add_argument("--mlow",dest="mlow", action="store", default=None,type=float,help="lower m1 and m2 component masses. Overrides mass intervals deduced from input argument")

  opts = parser.parse_args()

  flow = opts.flow
  fup  = opts.fup
  f_in = opts.input
  mlow = opts.mlow

  print("flow = %f, fup = %f"%(flow,fup))

  print mlow

  raise ValueError("script does not provide point at fup... fix this for future runs")

  # this will give a conservative estimate for df since m1,m2 are minimized independently
  for f in f_in:
    if mlow is None:
      param_intervals = np.genfromtxt(f,names=True)
      print("file %s with param intervals like..."%f)
      print param_intervals
      Mcs  = param_intervals['chirpmass'][0:2]
      etas = param_intervals['eta'][0:2]
      Mcs = np.linspace(Mcs.min(), Mcs.max(),150)
      etas = np.linspace(etas.min(), etas.max(),150)
      m1_min = 10000000.
      m2_min = 10000000.
      for Mc in Mcs:
        for eta in etas:
          m1, m2 = Mcnu_to_m1m2(Mc, eta)
          if m1 > 1. and m2 > 1.:
            #print("m1 = %f, m2 = %f"%(m1,m2))
            if m1 < m1_min:
              m1_min = m1
            if m2 < m2_min:
              m2_min = m2
    else:
      m1_min = mlow
      m2_min = mlow

    print("Lower m1 = %f, m2 = %f"%(m1_min,m2_min))

    nodes, weights = GravitationalWaveAdaptive(flow,fup,m1_min,m2_min)

    prefix = f.split('.')[0]
    np.savetxt(prefix+'_Nodes.dat', nodes)
    np.savetxt(prefix+'_Weights.dat', weights)



main()
