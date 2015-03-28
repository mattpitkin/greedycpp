""" collection of common python functions needed 
for analysis of greedy's output"""

import numpy as np

def load_ts(outdir):

    TS_real = np.loadtxt(outdir+'/TSpace_real.txt')
    TS_imag = np.loadtxt(outdir+'/TSpace_imag.txt')
    evaluations, quad_points = TS_real.shape
    TS = np.zeros((quad_points,evaluations), dtype=np.complex)
    TS[:,:].real = TS_real.transpose()
    TS[:,:].imag = TS_imag.transpose()

    return TS

def load_quad(outdir):
    quad_rule  = np.loadtxt(outdir+'/quad_rule.txt')
    weights    = quad_rule[:,1]
    nodes      = quad_rule[:,0]
    return weights,nodes

def load_info(outdir):

    ### load the basis and quadrature weights ###
    basis_real = np.loadtxt(outdir+'/Basis_real.txt')
    basis_imag = np.loadtxt(outdir+'/Basis_imag.txt')
    quad_rule  = np.loadtxt(outdir+'/quad_rule.txt')
    weights    = quad_rule[:,1]
    nodes      = quad_rule[:,0]

    ### remove the zeros --- final XXX rows could be identically zero ###
    basis_real = basis_real[np.where(np.any(basis_real != 0, axis=1))]
    basis_imag = basis_imag[np.where(np.any(basis_imag != 0, axis=1))]

    bases, quad_points = basis_real.shape
    B = np.zeros((quad_points, bases), dtype=np.complex)

    B[:,:].real = basis_real.transpose()
    B[:,:].imag = basis_imag.transpose()

    return B, weights, nodes
