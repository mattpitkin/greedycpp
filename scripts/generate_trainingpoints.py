#!/usr/bin/env python
"""
Created on Tue June 23 3:30:48 2015

@author: Jeroen Meidam
"""


"""
Generate a trainingpoints file for a veriety and combination of parameter options
"""

usage="""  generate_trainingpoints.py --input input.txt [options]
  This script is based on the generate_param_samples.py, but provides a bit more
  freedom in choice of parameters and choice of how to vary certain parameters.

  It reeds in a file where each column represents min,max,dimension and sampling type respectively.
  An example file could be one where you vary deterministically uniformly (\"linear\")
  over neutron star masses (col1&2) and randomly uniformly over component spin magnitudes along the z-axis (col5&8):
    1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0
    3.0 3.0 0.0 0.0 1.0 0.0 0.0 1.0
    50 50 1 1 50 1 1 50
    det_lin det_lin none none rand_lin none none rand_lin
  This will create a trainingspace of dimension 50x50x50x50. The parameters with dim=1
  will be fixed to the min value

  Extending an existing basis with an additional parameter using the --extend option
  The extension file provided with --extend is treated in a special way. The
  File could for example contain values for m1 and m2 which you would like to extend
  with a spin component sampled uniformly between -1 and 1 with 50 samples. The basis file
  could look like
    #m1 m2
    1.2 1.4
    1.3 1.3
    2.0 1.8
  The input file would be
    0. 0. -1.
    0. 0.  1.
    1  1  50
    file_0 file_1 rand_lin
  The resulting file would then have a dimension of 3*50, regardless of how many columns there
  were in the basis file; to each spin value the complete existing basis is added.
  The 0 and 1 index indicate which column in the basis file to use

  The --read-from-vectors options allows you to provide a predefined sampled set of points for
  a certain parameter. In the first case, in stead of sampling mass1 and mass2 uniformly you could use
  predefined vectors of e.g. length 50:
    1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0
    3.0 3.0 0.0 0.0 1.0 0.0 0.0 1.0
    50 50 1 1 50 1 1 50
    vector_0 vector_1 none none rand_lin none none rand_lin
  Now the resulting file would again be 50x50x50x50 in length but instead of automatically generating
  a vector for mass1 and mass2, the files 0 and 1 were used. The indices refer to the files provided
  in a space separated list with e.g. --read-from-vectors file0.txt file1.txt.
  When specifying a vector file, min, max and dim values are irrelevant and ignored.

  It is also possible to estimate the memory requirements for the trainingspace
  if you set --estimate-memory. In that case the following three quantities
  need to be provided: --max-rb, --quad-size and --num-workers.
"""

import numpy as np
import math, random, time
from optparse import OptionParser

### memory calculation function is copied from the memory module (this has no print statements)
def estimate_memory(ts_size,param_dim,max_rb,num_workers,quad_size):

    sizeof_double = 8.0e-9 # double in GB
    memory_worker = sizeof_double*( ts_size*param_dim + num_workers + 2*max_rb*ts_size/num_workers + 2*quad_size*ts_size/num_workers )
    memory_master = 2*sizeof_double*( num_workers + quad_size*max_rb + max_rb*max_rb) + sizeof_double*ts_size*param_dim

    return memory_worker*num_workers + memory_master*5.0

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

def CaclulateIndices_FileInput(iteration,params_num,params,file_length):
  """
  For use in generate_trainingpoints
  Each parameter has a vector of a certain length and the current iteration
  needs to pick out the correct number from each of these vectors.
  This function returns the indices of the parameter vectors that correspond
  to the current iteration.
  For example, take 1 parameters, with vector [1,2] and a file that
  has several columns of length 3. This will
  result in the following list:
  0 0 0 0 1
  1 0 0 0 2
  2 1 1 1 1
  3 1 1 1 2
  4 2 2 2 1
  5 2 2 2 2
  At iteration 4, indices = [2,2,2,0] and at iteration 1 indices = [0,0,0,1]
  """
  NP = len(params)

  #check which parameters have to come from a file
  files = []
  dims = [file_length]
  for i in range(NP):
    if params[i].sampling_mode == "file":
      files.append(i)
    if i > 0:
      dims.append(params[i].dim)

  def fn_k(n,k,dims):
    prod1 = np.prod(dims[:k+1])
    prod2 = 1.0
    if k > 0:
      prod2 = np.prod(dims[:k])
    return (n%prod1)/prod2

  indices = []
  for k in range(len(dims)):
    if k in files:
      fnk = iteration%file_length
    else:
      fnk = fn_k(iteration,k,dims)
    indices.append(fnk)

  return indices

def index_test():
  Nf = 6
  N2 = 4
  N4 = 2
  Ntotal = Nf*N2*N4
  for n in range(Ntotal):
    print n%Nf,n%Nf,n%N2,0,n%(Nf*N2)

def parse_samplingtype_string(sampling):
  sampling_mode = "deterministic"
  sampling_type = "linear"

  modes = {"rand":"random","det":"deterministic"}
  types = {"lin":"linear","ln":"ln","log10":"log10","power":"power"}

  if sampling.strip() != "none":
    spl = sampling.strip().split('_')
    print spl
    print len(spl)
    if len(spl) not in [2,3]:
      raise Exception("Could not parse string \"%s\"."%sampling)
    else:
      sampling_mode = spl[0]
      sampling_type = spl[1]
      if modes.has_key(spl[0]):
        sampling_mode = modes[spl[0]]
      elif sampling_mode == "file":
        pass
      elif sampling_mode == "vector":
        pass
      else:
        raise Exception("Undefined sampling mode \"%s\" (use \"det\" or \"rand\")."%spl[0])

      if types.has_key(spl[1]):
        sampling_type = types[spl[1]]
      elif sampling_mode == "file":
        sampling_type = spl[1]
      elif sampling_mode == "vector":
        sampling_type = spl[1]
      else:
        raise Exception("Undefined sampling type \"%s\" (use \"lin\", \"power\", \"ln\" or \"log10\")."%spl[1])

      if sampling_type in ['power', 'log10', 'ln']:
        log_sampling_factor = float(spl[2])
        print("using sample factor %f"%log_sampling_factor)
      else:
        log_sampling_factor = None
        

  return [sampling_mode, sampling_type, log_sampling_factor]


class parameter(object):
  """ Parameter object storing information for each parameter """
  def __init__(self,parameter_index,min,max,dim,sampling,data_to_extend=[],vector_files=[]):

    self.min = min # minimum value
    self.max = max # maximum value
    self.dim = dim # number of samples from min to max
    smpl = parse_samplingtype_string(sampling)
    self.sampling_mode = smpl[0]
    self.sampling_type = smpl[1]
    self.scaling_factor = smpl[2] # for logarithmic sampling
    self.index = parameter_index

    if dim < 1:
      raise Exception("Parameter dimensions must be at least 1")

    if self.sampling_mode == "file":
      if len(data_to_extend)==0:
        raise Exception("No data file provided. Data is empty.")
      self.dim = 1
      print "p%d: data from file, column %s (length=%d)"%(self.index,self.sampling_type,len(data_to_extend[:,0]))
      self.vector = data_to_extend[:,int(self.sampling_type)]
    elif self.sampling_mode == "vector":
      vfile = vector_files[int(self.sampling_type)]
      vector_from_file=np.genfromtxt(vfile)
      self.dim = len(vector_from_file)
      print "p%d: using vector from file \"%s\" (length=%d)"%(self.index,vfile,self.dim)
      self.vector = vector_from_file
    else:
      if dim > 1:
        print "p%d: sample_mode %s, sample_type %s (min=%.2f,max=%.2f,dim=%d)"%(self.index,self.sampling_mode,self.sampling_type,self.min,self.max,self.dim)
        self.vector = self.generate_vector()
      else:
        print "p%d: constant (value=%.2f)"%(self.index,self.min)
        self.vector = np.array([self.min])

  def generate_vector(self):

    if( self.sampling_mode is "deterministic"):
      if(self.sampling_type=="linear"):
        v = np.linspace(self.min,self.max,self.dim,endpoint=True)
      elif(self.sampling_type=="ln"):
        v_unscaled = np.logspace(math.log(self.min), math.log(self.scaling_factor*self.max), self.dim, endpoint=True, base=math.exp(1))
        v=self._rescale_log(v_unscaled, self.min, self.max, self.scaling_factor)
      elif(self.sampling_type=="log10"):
        v_unscaled = np.logspace(math.log10(self.min), math.log10(self.scaling_factor*self.max), self.dim, endpoint=True, base=10.0)
        v=self._rescale_log(v_unscaled, self.min, self.max, self.scaling_factor)
      elif(self.sampling_type=="power"):
        v=self._power_points(self.scaling_factor, self.min, self.max, self.dim)

    elif ( self.sampling_mode is "random"):
      v = np.array([(self.max - self.min) * random.random() + self.min for i in range(self.dim)])

    else:
      raise Exception("sampling mode unknown")

    return v

  def _rescale_log(self, x, p_min, p_max, scaling_factor):
    "rescale the log range to be in [p_min,p_max] if upper log-range has been pushed out"
    a,b0,b = p_min,p_max,p_max*scaling_factor
    c0,c1=a*(-b + b0)/(a - b),(a - b0)/(a - b)
    return c0+c1*x

  def _power_points(self, p, var_min, var_max, N):
    var_p = np.linspace(var_min**p, var_max**p, N)
    return var_p**(1./p)

def generate_trainingpoints(filename,params_low,params_high,params_num,params_type,file_to_extend="",vector_files=[], return_ts=False):

  Nparams = len(params_low)
  if (Nparams != len(params_high) or Nparams != len(params_num) or Nparams != len(params_type)):
    raise Exception("dimensions of params_low, params_high, params_num, params_type do not agree (%d,%d,%d,%d)"%(Nparams,len(params_high),len(params_num),len(params_type)))

  #determine the training space size N and create parameter objects
  Nfile=1
  if file_to_extend=="":
    params = [parameter(0,params_low[0],params_high[0],params_num[0],params_type[0],vector_files=vector_files)]
    N = params[0].dim
    for i in range(1,Nparams):
      params.append( parameter(i,params_low[i],params_high[i],params_num[i],params_type[i],vector_files=vector_files) )
      N *= params[i].dim
  else:
    data = np.genfromtxt(file_to_extend,dtype=float)
    Nfile = len(data[:,0])
    #N = Nfile #TS size is at least the length of the file to extend
    params = [parameter(0,params_low[0],params_high[0],params_num[0],params_type[0],data_to_extend=data,vector_files=vector_files)]
    N = params[0].dim
    for i in range(1,Nparams):
      params.append( parameter(i,params_low[i],params_high[i],params_num[i],params_type[i],data_to_extend=data,vector_files=vector_files) )
      N *= params[i].dim

  tic = time.time()
  Ntotal = N*Nfile
  print "Generating %d trainingpoints..."%(Ntotal)
  if file_to_extend=="":
    with open(filename,"w") as f:
      # much faster to precompute all pv_indicies
      #pv_indices = CaclulateIndices(0,params_num,len(params))
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
  else:
    with open(filename,"w") as f:
      for n in range(Ntotal):#for n in range(file_length*N):
        param_vector_indices = CaclulateIndices_FileInput(n,params_num,params,Nfile)
        string = ""
        for i in range(Nparams):
          string+=str(params[i].vector[param_vector_indices[i]])
          if i < Nparams-1:
            string+=" "
          if i == Nparams-1 and n < Ntotal-1:
            string+="\n"

        f.write(string)

  toc = time.time()
  print("Seconds to sample: %f"%(toc-tic))

  if return_ts: 
    ts = np.loadtxt(filename)
  else:
    ts = None

  return Ntotal, ts



def generate_sampling(filename,param_list=None,total_picks=100,seed=None,file_to_extend="",vector_files=[],return_ts=False):
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

  tic = time.time()
  TS_size,TS = generate_trainingpoints(filename,params_low,params_high,params_num,params_type,file_to_extend=file_to_extend,vector_files=vector_files,return_ts=return_ts)
  toc = time.time()
  print("Seconds to finish generate_trainingpoints: %f"%(toc-tic))
 

  return TS_size, TS

def vararg_callback(option, opt_str, value, parser):
    assert value is None
    value = []
    def floatable(str):
        try:
            float(str)
            return True
        except ValueError:
            return False
    for arg in parser.rargs:
        # stop on --foo like options
        if arg[:2] == "--" and len(arg) > 2:
            break
        # stop on -a, but not on -3 or -3.0
        if arg[:1] == "-" and len(arg) > 1 and not floatable(arg):
            break
        value.append(arg)
    del parser.rargs[:len(value)]
    setattr(parser.values, option.dest, value)


def main():

  parser=OptionParser(usage)
  parser.add_option("-i","--input",dest="input", action="store", type="string",help="input file containing sampling ranges, modes and types")
  parser.add_option("-o","--out",dest="fout", action="store", default="./trainingpoints.txt", type="string",help="output file to store the trainingpoints")
  parser.add_option("-s","--seed",dest="seed", action="store", default=1234, type=int,help="seed for random generation")
  parser.add_option("-m","--estimate-memory",default=False,dest="estimate_memory",action="store_true",help="Estimate memory requirements for trainingspace")
  parser.add_option("-M","--estimate-memory-only",default=False,dest="mem_only",action="store_true",help="Estimate only memory")
  parser.add_option("-q","--quad-size",dest="quad_size", action="store", type=int,help="length of quads")
  parser.add_option("-r","--max-rb",dest="max_rb", action="store", type=int, help="max reduced basis size that is set in the config file")
  parser.add_option("-t","--size-ts",dest="ts_size", action="store", type=int, help="Size of training set, required if only calculating memory")
  parser.add_option("-w","--num-workers",dest="num_workers", default=1, action="store", type=int,help="number of workers used")
  parser.add_option("-e","--extend-file",dest="extend_file", action="store", type="string",help="extend the input file with parameters from input")
  parser.add_option("-v","--read-from-vectors",dest="vectors",action="callback",callback=vararg_callback,help="If a column's type is \"vector\" use this argument to pass the vector files (space separated)")
  (opts,args) = parser.parse_args()


  if not opts.mem_only:
    if not opts.input:
      print "Error: Need input file containing parameter information. Specify with --input.\n"
      exit(usage)
    else:
      f_in = opts.input

  random.seed(opts.seed)

  f_out   = opts.fout

  f_extend=""
  if opts.extend_file:
    f_extend = opts.extend_file

  vectors=[]
  if opts.vectors:
    vectors = opts.vectors

  if opts.mem_only:
    if not opts.max_rb:
      exit("Error in estimating memory: must provide size of training set if --estimate-memory-only is set")
    ts_size = opts.ts_size
  else:
    ts_size, ts = generate_sampling(f_out,param_list=f_in,file_to_extend=f_extend,vector_files=vectors)
    print "Trainingpoints written to \"%s\""%f_out

  ### Perform memory estimation if requested ###
  if (opts.estimate_memory) or (opts.mem_only):
    if not opts.max_rb:
      exit("Error in estimating memory: Must provide --max-rb and --quad-size for memory estimate")
    if not opts.quad_size:
      exit("Error in estimating memory: Must provide --max-rb and --quad-size for memory estimate")

    mem = estimate_memory(ts_size,1,opts.max_rb,opts.num_workers,opts.quad_size)
    print 'Run requires         ~%.2f GB' % mem
    print '  10 pc fudge factor: %.2f GB' % (mem*1.1)
    print '  20 pc fudge factor: %.2f GB' % (mem*1.2)


if __name__ == "__main__":
	# START THE MAIN FUNCTION IF RUN AS A SCRIPT. OTHERWISE, JUST LOADS THE CLASS
	# AND FUNCTION DEFINITIONS

        tic = time.time()
	main()
        toc = time.time()    
        print("Seconds to run sampler: %f"%(toc-tic))
        exit()

