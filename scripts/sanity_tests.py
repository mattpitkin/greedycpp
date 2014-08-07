"""
Module to do sanity tests of output directory...

1) basis orthogonality test
2) 


Can run at command line from an application directory, e.g.:
  $ python sanity_tests.py OUTDIR1
"""

import os, re, sys
import numpy as np


def ortho_test(outdir):

    ### load the basis and quadrature weights ###
    basis_real = np.loadtxt(outdir+'/Basis_real.txt')
    basis_imag = np.loadtxt(outdir+'/Basis_imag.txt')
    weights    = np.loadtxt(outdir+'/quad_weights.txt')

    ### remove the zeros --- final XXX rows are identically zero ###
    basis_real = basis_real[np.where(np.any(basis_real != 0, axis=1))]
    basis_imag = basis_imag[np.where(np.any(basis_imag != 0, axis=1))]

    bases, quad_points = basis_real.shape
    B = np.zeros((quad_points, bases), dtype=np.complex)

    B[:,:].real = basis_real.transpose()
    B[:,:].imag = basis_imag.transpose()

    result = weights[0] * np.dot(B.T.conj(),B) # NOTE: ASSUMES ALL WEIGHTS ARE THE SAME (TODO)

    err = result - np.eye(result.shape[1])
    print "othogonality error %1.15e" % np.max(np.abs(err))

if __name__=="__main__":

    try:
        outdir = sys.argv[1]
    except:
        raise Exception("Must specify an output directory")

    ortho_test(outdir)

