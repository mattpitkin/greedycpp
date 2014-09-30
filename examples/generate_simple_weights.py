'''from a file of quadrature points, generate set of weights for Reimann sum'''

import numpy as np
import sys

def GenerateWeights(filename):

    quad_pts = np.loadtxt(filename)
    weights  = np.diff(quad_pts)
    weights  = np.append(weights,0)

    fp = open('MyWeights.txt','w')
    for ii in range(np.size(weights)):
        fp.write('%1.15e\n' % weights[ii])
    fp.close() 


if __name__=="__main__":

    print sys.argv[1]

    try:
        filename = sys.argv[1]
    except:
        raise Exception("Must specify a file of quadrature points")

    GenerateWeights(filename)
