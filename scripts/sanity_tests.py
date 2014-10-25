"""
Module to do sanity tests of output directory...

1) basis orthogonality test
2) basis accuracy test (training space file needs to be present)


Can run at command line from an application directory, e.g.:
  $ python sanity_tests.py OUTDIR1
"""

import os, re, sys
import numpy as np

def load_info(outdir):

    ### load the basis and quadrature weights ###
    basis_real = np.loadtxt(outdir+'/Basis_real.txt')
    basis_imag = np.loadtxt(outdir+'/Basis_imag.txt')
    quad_rule  = np.loadtxt(outdir+'/quad_rule.txt')
    weights    = quad_rule[:,1]
    nodes      = quad_rule[:,0]

    ### remove the zeros --- final XXX rows are identically zero ###
    basis_real = basis_real[np.where(np.any(basis_real != 0, axis=1))]
    basis_imag = basis_imag[np.where(np.any(basis_imag != 0, axis=1))]

    bases, quad_points = basis_real.shape
    B = np.zeros((quad_points, bases), dtype=np.complex)

    B[:,:].real = basis_real.transpose()
    B[:,:].imag = basis_imag.transpose()

    return B, weights, nodes

def load_ts(outdir):

    TS_real = np.loadtxt(outdir+'/TSpace_real.txt')
    TS_imag = np.loadtxt(outdir+'/TSpace_imag.txt')
    evaluations, quad_points = TS_real.shape
    TS = np.zeros((quad_points,evaluations), dtype=np.complex)
    TS[:,:].real = TS_real.transpose()
    TS[:,:].imag = TS_imag.transpose()

    return TS

def ortho_test(outdir):

    B, weights, nodes = load_info(outdir)

    result = np.dot(weights * B.T.conj(),B) 

    err = result - np.eye(result.shape[1])
    print "othogonality error %1.15e" % np.max(np.abs(err))

def basis_accuracy_test(outdir):

    B, weights, nodes = load_info(outdir)

    TS = load_ts(outdir)
    quad_points, evaluations = TS.real.shape 

    err_app = np.zeros(evaluations)

    # TODO: shouldn't use for-loop here
    for ii in range(evaluations):
        proj_coeff = np.dot(weights * B.T.conj(),TS[:,ii])
        err_h = TS[:,ii] - np.dot(B,proj_coeff) # not conj is good! we just want the sum (no weights either)
        err_app[ii] = np.sqrt( np.dot( weights * err_h,err_h.conj()) ).real

    import matplotlib.pyplot as plt
    plt.semilogy(range(evaluations),err_app)
    plt.show()

    ## to show other info ##
    #import matplotlib.pyplot as plt
    #plt.plot(nodes,TS[:,1])
    #plt.show()
    #plt.plot(nodes,err_h)
    #plt.show()


if __name__=="__main__":

    try:
        outdir = sys.argv[1]
    except:
        raise Exception("Must specify an output directory")

    ortho_test(outdir)
    basis_accuracy_test(outdir)

