''' generate a mutlidomain quadrature rule for input into the greedy alg '''

import numpy as np
from scipy.special.orthogonal import p_roots


### setup quadrature rule ###
pts = 10000
#dom_ints = [10.0,30.0,1024.0]
dom_ints = [40,65,100,250,400,700,1099]

### make these global to avoid recomputing roots ###
p_roots_x, p_roots_w = p_roots(pts)

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

quad_nodes,quad_weights = MultiDomain(pts,dom_ints)

fp = open('MyQuadPoints.txt','w')
fw = open('MyQuadWeights.txt','w')
for ii in range(len(quad_nodes)):
    fp.write('%1.15e\n' % quad_nodes[ii])
    fw.write('%1.15e\n' % quad_weights[ii])

fp.close()
fw.close()


