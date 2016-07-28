#!/usr/bin/env python

''' generate a mutlidomain quadrature rule for use with the greedy algorithm. '''

import numpy as np
from scipy.special.orthogonal import p_roots


points_fname  = 'MyQuadPoints.txt'
weights_fname = 'MyQuadWeights.txt'


#quadrature_type = 'GL'
#quadrature_type = 'MultiDomain'
quadrature_type = 'GW'

### setup quadrature rule ###
#pts = 10000
##dom_ints = [10.0,30.0,1024.0]
#dom_ints = [40,65,100,250,400,700,1099]

### make these global to avoid recomputing roots ###
#p_roots_x, p_roots_w = p_roots(pts)

# for use with the routine GravitationalWaveAdaptive
fmin = 20
fmax = 4096
dflim = 5
m1_min = 1.2
m2_min = 1.2

def SingleDomGaussLeg(quad_points,a,b):
    '''a,b domain endpoints. quad_points is number of quadrature points'''

    #x,w = p_roots(quad_points)
    x = p_roots_x
    w = p_roots_w
    x = x.real
    y = (b-a)*(x+1.0)/2.0 + a
    w = (b-a)/2.0*w
    return y,w

def MultiDomain(quad_points,dom_intervals):
    '''[x0,x1],[x1,x2]... are the domain intervals'''

    for i in range(1,len(dom_intervals)):
      if ( dom_intervals[i] <= dom_intervals[i-1]): 
          raise Exception,"boundaries of intervals must be monotonic."

    x = np.zeros((0,0),float)
    w = np.zeros((0,0),float)
    for i in range(len(dom_intervals)-1):
        x_dom,w_dom = SingleDomGaussLeg(quad_points,dom_intervals[i],dom_intervals[i+1])
        x = np.append(x,x_dom)
        w = np.append(w,w_dom)

    return x,w

def GravitationalWaveAdaptive(fmin,fmax,dflim,m1_min,m2_min):
  ''' Simple algorithm to create an array of frequencies and deltaF's
      based on gravitational wave signal duration.

      Written by Michael Purrer.

  INPUT
  =====
  fmin   --- starting frequency 
  fmax   --- ending frequency
  dflim  --- maximum allowed df
  m1_min --- smallest component mass of body 1
  m2_min --- smallest component mass of body 2'''

  import lal
  import lalsimulation as LS

  df_array = []
  f_array = []
  f = fmin

  while f < fmax:
    try:
      df = 0.2 / LS.SimIMRSEOBNRv2ROMDoubleSpinHITimeOfFrequency(f, m1_min*lal.MSUN_SI, m2_min*lal.MSUN_SI, 0, 0)
    except Exception,e:
      print str(e)
      df = dflim # At very high frequencies the above call can fail, but we know that the frequency spacing would be huge there
    #print f, df
    if abs(df) > dflim: df = dflim
    f_array.append(f)
    df_array.append(df)
    f += df

  print len(f_array)

  return f_array, df_array



if quadrature_type == 'MultiDomain':
  quad_nodes, quad_weights = MultiDomain(pts,dom_ints)
elif quadrature_type == 'GW':
  quad_nodes, quad_weights = GravitationalWaveAdaptive(fmin,fmax,dflim,m1_min,m2_min)
elif quadrature_type == 'GL':
  #SingleDomGaussLeg(quad_points,a,b)
  pass
else:
  raise ValueError

np.savetxt(points_fname, quad_nodes)
np.savetxt(weights_fname, quad_weights)


#fp = open(points_fname,'w')
#fw = open(weights_fname,'w')
#for ii in range(len(quad_nodes)):
  #fp.write('%1.15e\n' % quad_nodes[ii])
  #fw.write('%1.15e\n' % quad_weights[ii])
#fp.close()
#fw.close()


