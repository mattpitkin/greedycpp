#!/usr/bin/env python
"""
Created on Tue June 23 3:30:48 2015

@author: Jeroen Meidam
"""


"""
Generate a trainingpoints file for a veriety and combination of parameter options
"""

usage="""  generate_trainingpoints.py input.txt [options]
  This script is based on the generate_param_samples.py, but provides a bit more
  freedom in choice of parameters and choice of how to vary certain parameters.

  It reeds in a file where each row represents min,max,dimension and sampling type respectively.
  An example file could be one where you vary deterministically uniformly (\"linear\")
  over neutron star masses (col1&2) and randomly uniformly over component spin magnitudes along the z-axis (col5&8):
    1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0
    3.0 3.0 0.0 0.0 1.0 0.0 0.0 1.0
    50 50 1 1 50 1 1 50
    det_lin det_lin none none rand_lin none none rand_lin
  This will create a trainingspace of dimension 50x50x50x50. The parameters with dim=1
  will be fixed to the min value

  It is also possible to estimate the memory requirements for the trainingspace
  if you set --estimate-memory. In that case the following three quantities
  need to be provided: --max-rb, --quad-size and --num-workers.
"""

import numpy as np
import math
import random
from optparse import OptionParser

### import memory.py from scripts directory
import memory

def CaclulateIndices(iteration,params_num,n_params):
  """
  For use in generate_trainingpoints
  Each parameter has a vector of a certain length and the current iteration
  needs to pick out the correct number from each of these vectors.
  This function returns the indices of the parameter vectors that correspond
  to the current iteration.
  For example, take 2 parameters, with vectors [1,2,3] and [1,2]. This will
  result in the following list:
  0 1 1
  1 1 2
  2 2 1
  3 2 2
  4 3 1
  5 3 2
  At iteration 4, indices = [2,0] and at iteration 1 indices = [0,1]
  """
  NP = n_params
  dims = params_num
  i = range(NP-1)
  i.reverse()
  indices=np.zeros(NP,dtype=int)
  indices[NP-1] = iteration
  for ii in i:
    indices[ii] = ( indices[ii+1] - indices[ii+1]%dims[ii+1] )/dims[ii+1]
  for ii in range(NP):
    indices[ii] = indices[ii]%dims[ii]

  return indices

def parse_samplingtype_string(sampling):
  sampling_mode = "deterministic"
  sampling_type = "linear"

  modes = {"rand":"random","det":"deterministic"}
  types = {"lin":"linear","ln":"ln","log10":"log10"}

  if sampling.strip() != "none":
    spl = sampling.strip().split('_')
    if len(spl)!=2:
      raise Exception("Could not parse string \"%s\"."%sampling)
    else:
      sampling_mode = spl[0]
      sampling_type = spl[1]
      if modes.has_key(spl[0]):
        sampling_mode = modes[spl[0]]
      else:
        raise Exception("Undefined sampling mode \"%s\" (use \"det\" or \"rand\")."%spl[0])
      if types.has_key(spl[1]):
        sampling_type = types[spl[1]]
      else:
        raise Exception("Undefined sampling type \"%s\" (use \"lin\", \"ln\" or \"log10\")."%spl[1])

  return [sampling_mode, sampling_type]


class parameter(object):
  """ Parameter object storing information for each parameter """
  def __init__(self,parameter_index,min,max,dim,sampling,scaling_factor = 20.0):

    self.min = min
    self.max = max
    self.dim = dim
    self.scaling_factor = scaling_factor
    smpl = parse_samplingtype_string(sampling)
    self.sampling_mode = smpl[0]
    self.sampling_type = smpl[1]
    self.index = parameter_index

    if dim < 1:
      raise Exception("Parameter dimensions must be at least 1")

    if dim > 1:
      print "p%d: sample_mode %s, sample_type %s (min=%.2f,max=%.2f,dim=%d)"%(self.index,self.sampling_mode,self.sampling_type,self.min,self.max,self.dim)
      self.vector = self.generate_vector()
    else:
      print "p%d: constant (value=%.2f)"%(self.index,self.min)
      self.vector = np.array([self.min])

  def rescale_log(x,p_min,p_max,scaling_factor):
    "rescale the log range to be in [p_min,p_max] if upper log-range has been pushed out"
    a,b0,b = p_min,p_max,p_max*scaling_factor
    c0,c1=a*(-b + b0)/(a - b),(a - b0)/(a - b)
    return c0+c1*x

  def generate_vector(self):

    if( self.sampling_mode is "deterministic"):
      if(self.sampling_type=="linear"):
        v = np.linspace(self.min,self.max,self.dim,endpoint=True)
      if(self.sampling_type=="ln"):
        v_unscaled = np.logspace(math.log(self.min), math.log(self.scaling_factor*self.max), self.dim, endpoint=True, base=math.exp(1))
        v=self.rescale_log(v_unscaled, self.min, self.max, self.scaling_factor)
      if(self.sampling_type=="log10"):
        v_unscaled = np.logspace(math.log10(self.min), math.log10(self.scaling_factor*self.max), self.dim, endpoint=True, base=10.0)
        v=self.rescale_log(v_unscaled, self.min, self.max, self.scaling_factor)

    elif ( self.sampling_mode is "random"):
      v = np.array([(self.max - self.min) * random.random() + self.min for i in range(self.dim)])

    else:
      raise Exception("sampling mode unknown")

    return v

def generate_trainingpoints(filename,params_low,params_high,params_num,params_type):

  Nparams = len(params_low)
  if (Nparams != len(params_high) or Nparams != len(params_num) or Nparams != len(params_type)):
    raise Exception("dimensions of params_low, params_high, params_num, params_type do not agree (%d,%d,%d,%d)"%(Nparams,len(params_high),len(params_num),len(params_type)))

  #determine the training space size N and create parameter objects
  N = params_num[0]
  params = [parameter(0,params_low[0],params_high[0],params_num[0],params_type[0])]
  for i in range(1,Nparams):
    N *= params_num[i]
    params.append( parameter(i,params_low[i],params_high[i],params_num[i],params_type[i]) )


  print "Generating %d trainingpoints..."%N
  with open(filename,"w") as f:
    for n in range(N):
      param_vector_indices = CaclulateIndices(n,params_num,len(params))
      string = ""
      for i in range(Nparams):
        string+=str(params[i].vector[param_vector_indices[i]])
        if i < Nparams-1:
          string+=" "
        if i == Nparams-1 and n < N-1:
          string+="\n"

      f.write(string)

  return N



def generate_sampling(filename,param_list=None,total_picks=100,seed=None):
  if param_list is None: # manual setup
    params_low  = np.array([1.0,1.0])  # lower interval of each parameter
    params_high = np.array([3.0,3.0])  # upper interval of each parameter
    params_num = np.array([50,50])  # vector length of each parameter
    params_type = np.array(["det_lin","det_lin"])
  else:
    try: # attempt to get the max/min from a file
      param_info = np.genfromtxt(param_list,dtype=str)
      params_low  = [float(p) for p in param_info[0]]
      params_high = [float(p) for p in param_info[1]]
      params_num = [int(p) for p in param_info[2]]
      params_type = param_info[3]
    except: # parameter ranges have been passed as lists
      param_info = param_list
      params_low  = param_info[0]
      params_high = param_info[1]
      params_num = param_info[2]
      params_type = param_info[3]

  TS_size = generate_trainingpoints(filename,params_low,params_high,params_num,params_type)

  return TS_size


def main():

  parser=OptionParser(usage)
  parser.add_option("-o","--out",dest="fout", action="store", default="./trainingpoints.txt", type="string",help="output file to store the trainingpoints")
  parser.add_option("-s","--seed",dest="seed", action="store", default=1234, type=int,help="seed for random generation")
  parser.add_option("-m","--estimate-memory",default=False,dest="estimate_memory",action="store_true",help="Estimate memory requirements for trainingspace")
  parser.add_option("-q","--quad-size",dest="quad_size", action="store", type=int,help="length of quads")
  parser.add_option("-r","--max-rb",dest="max_rb", action="store", type=int, help="max reduced basis size that is set in the config file")
  parser.add_option("-w","--num-workers",dest="num_workers", action="store", type=int,help="number of workers used")
  (opts,args) = parser.parse_args()

  if len(args) != 1:
    print "Error: Need 1 argument; input file containing parameter information\n"
    exit(usage)

  random.seed(opts.seed)

  f_in    = args[0]
  f_out   = opts.fout

  ts_size = generate_sampling(f_out,f_in)
  print "Trainingpoints written to \"%s\""%f_out

  ### Perform memory estimation if requested ###
  if (opts.estimate_memory):
    if not opts.max_rb:
      exit("Error in estimating memory: Must provide --max-rb,--num-workers and --quad-size for memory estimate")
    if not opts.num_workers:
      exit("Error in estimating memory: Must provide --max-rb,--num-workers and --quad-size for memory estimate")
    if not opts.quad_size:
      exit("Error in estimating memory: Must provide --max-rb,--num-workers and --quad-size for memory estimate")

    memory.estimate_memory(ts_size,1,opts.max_rb,opts.num_workers,opts.quad_size)
#    num_workers = opts.num_workers
#    max_rb = opts.max_rb
#    quad_size = opts.quad_size
#    sizeof_double = 8.0e-9 # double in GB
#    memory_worker = sizeof_double*( ts_size + num_workers + 2*max_rb*ts_size/num_workers + 2*quad_size*ts_size/num_workers )
#    memory_master = 2*sizeof_double*( num_workers + quad_size*max_rb + max_rb*max_rb) + sizeof_double*ts_size
#    print 'All worker processor require at least %1.2f GB (in total)' % (memory_worker*num_workers)
#    print 'Master (orthogonalization) processor requires at least %1.2f GB' % memory_master



if __name__ == "__main__":
	# START THE MAIN FUNCTION IF RUN AS A SCRIPT. OTHERWISE, JUST LOADS THE CLASS
	# AND FUNCTION DEFINITIONS
	exit(main())
