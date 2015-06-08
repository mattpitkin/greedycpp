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

def load_quad(outdir=None,xmin=None,xmax=None,quad_points=None):
  ''' if outdir is supplied, the quadrature rule is loaded.
      otherwise x=np.linspace(xmin,xmax,quad_points) is used'''

  if(outdir is not None):
    quad_rule  = np.loadtxt(outdir+'/quad_rule.txt')
    weights    = quad_rule[:,1]
    nodes      = quad_rule[:,0]
  else:
    nodes   = np.linspace(xmin,xmax,quad_points)
    weights = (nodes[1]-nodes[0])*np.ones(nodes.shape)
  return weights, nodes

def load_basis(outdir):

  ### load the basis and quadrature weights ###
  try: 
    basis_real = np.loadtxt(outdir+'/Basis_real.txt')
    basis_imag = np.loadtxt(outdir+'/Basis_imag.txt')

    ### remove the zeros --- final XXX rows could be identically zero ###
    basis_real = basis_real[np.where(np.any(basis_real != 0, axis=1))]
    basis_imag = basis_imag[np.where(np.any(basis_imag != 0, axis=1))]

    bases, quad_points = basis_real.shape
    B = np.zeros((quad_points, bases), dtype=np.complex)

    B[:,:].real = basis_real.transpose()
    B[:,:].imag = basis_imag.transpose()

    print 'Found basis in a text file'

    return B
  except IOError: 
    try:
      B = np.load(outdir+'/Basis.npy')
      B = B.transpose()
      print 'Found basis in a numpy file'
      return B
    except IOError:
      print 'No valid basis files could be found'

def load_eim(outdir):

  ### load the EIM point information ###
  eim_indx  = np.loadtxt(outdir+'/EIM_indices.txt').astype(int)
  eim_nodes = np.loadtxt(outdir+'/EIM_nodes.txt')

  ### load the EIM inverse matrix ###
  try:

    invV = np.zeros((len(eim_indx),len(eim_indx)), dtype=np.complex)

    invV_real      = np.loadtxt(outdir+'/invV_real.txt')
    invV_imag      = np.loadtxt(outdir+'/invV_imag.txt')
    invV[:,:].real = invV_real
    invV[:,:].imag = invV_imag

    print 'Found eim invV in a text file'

  except IOError:
    try:
      invV = np.load(outdir+'/invV.npy')
      print 'Found eim invV in a numpy file'

    except IOError:
      print 'No valid basis files could be found'
      raise IOError

  return invV, eim_indx, eim_nodes

def load_info(outdir,xmin=None,xmax=None,quad_points=None):
  ''' if outdir is supplied, the quadrature rule is loaded.
      otherwise x=np.linspace(xmin,xmax,quad_points) is used'''

  B = load_basis(outdir)
  if(xmax is not None):
    weights, nodes = load_quad(None,xmin,xmax,quad_points)
  else:
    weights, nodes = load_quad(outdir)

  return B, weights, nodes
