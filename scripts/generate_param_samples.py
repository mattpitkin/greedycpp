''' Generate parameter sampling, write these values to file. The resulting
    text file can be read as input to greedy.cpp through proper specification
    in the configuration file.'''

import numpy as np
import math
import sys, os
import timeit

# determine max/min of each parameter range (e.g. by running
# on a file of greedy points)
def get_param_max_min(filename):

  params = np.loadtxt(filename)
  param_range = np.zeros([2,params.shape[1]])
  for i in range(params.shape[1]):
    param_range[0,i] = params[:,i].min()
    param_range[1,i] = params[:,i].max()
    print "param %i max %1.14e and min %1.14e"%(i,param_range[1,i],param_range[0,i])

  return param_range

def generate_sampling(filename,param_list=None,total_picks=100,\
                      sample_type="rand",seed=None,return_picks=False):
    """ filename     -- name of output file
        param_list   -- optional. If an input file, a min/max range
                        of parameters is deduced. Otherwise, a vector 
                        of max/min parameter values directly.
        total_picks  -- if random sampling, this many random picks
        sample_type  -- 'rand' or 'deterministc' 
        seed         -- if 'rand', seed the RandomState class"""

    ### setup parameter intervals ###
    ### ith parameter interval will be [param_low[i],param_high[i]]
    if param_list is None: # manual setup
      params_low  = np.array([1.0,1.0])  # lower interval of each parameter
      params_high = np.array([3.0,3.0])
      #params_low  = np.array([2.8,0.098765,-.7,-.7,0.0,0.0])  # lower interval of each parameter
      #params_high = np.array([20.0,0.25,.7,-0.046667,2*np.pi,2*np.pi])
    else:
      try: # attempt to get the max/min from a file
        param_range = get_param_max_min(param_list)
        params_low  = param_range[0,:]
        params_high = param_range[1,:]
      except: # parameter ranges have been passed
        param_range = param_list
        params_low  = param_range[0,:]
        params_high = param_range[1,:]

    ### setup for deterministic sampling ###
    param_sampling = "ln"
    scaling_factor = 20.0 # default should be 1. higher values densely sample lower range of parameter interval
    params_num  = [50,50] # deterministic: upper interval of each parameter



    if( sample_type is "deterministic"):
      print "deterministic sampling"

      ### from franks code ###
      def rescale_mass(x,p_min,p_max):
        "rescale the mass range to be in [m_min,m_max] if upper log-range has been pushed out"
        a,b0,b = p_min,p_max,p_max*scaling_factor
        c0,c1=a*(-b + b0)/(a - b),(a - b0)/(a - b)
        return c0+c1*x

      p_low  = params_low[0]
      p_high = params_high[0]
      p_num  = params_num[0]
      if(param_sampling=="linear"):
        p1 = np.linspace(p_low,p_high,p_num,endpoint=True)
      if(param_sampling=="ln"):
        mass_params_unscaled = np.logspace(math.log(p_low), math.log(scaling_factor*p_high), p_num, endpoint=True, base=math.exp(1))
        mass_params=rescale_mass(mass_params_unscaled, p_low, p_high)
        p1 = mass_params
      if(param_sampling=="log10"):
        mass_params_unscaled = np.logspace(math.log10(p_low), math.log10(scaling_factor*p_high), p_num, endpoint=True, base=10.0)
        mass_params=rescale_mass(mass_params_unscaled)
        p1 = mass_params

      ### this could be more explict use of log spacing ###
      #p1 = np.power(params_low[0]*(params_high[0]/params_low[0]),np.linspace(0,1,params_num[0]))
      #p2 = np.power(params_low[1]*(params_high[1]/params_low[1]),np.linspace(0,1,params_num[1]))

      p2 = p1

      ### sanity checks ###
      if p1.min()<p_low: raise Exception("samples too low!")
      if p1.max()>p_high: raise Exception("samples too high!")
      for i in range(1,len(p1)):
          if p1[i]-p1[i-1]<=0: raise Exception("samples not monotonic!")

      if(return_picks):
        raise ValueError('Not coded yet')

      ### output MyTS.txt here -- tensor product parameter grid ###
      fp = open(filename,'w')
      for ii in range(np.size(p1)):
        for jj in range(np.size(p2)):
          fp.write('%1.15e\t' % p1[ii])
          fp.write('%1.15e\n' % p2[jj])

      fp.close()

    elif(sample_type is "rand"):

      rs = np.random.RandomState(seed=seed)

      fp = open(filename,'w')
      parameter_dim = len(params_high)

      tic = timeit.default_timer()
      for ii in range(total_picks):

        p_jj = (params_high[:] - params_low[:]) * rs.random_sample((parameter_dim,)) + params_low[:]

        if np.mod(ii,100000)==0:
          print "sample number = ",ii
          print "ii/total_picks  = %1.14f"%(float(ii)/float(total_picks))

        for jj in range(parameter_dim):
          if(jj == parameter_dim-1):
            fp.write('%1.15e\n' % p_jj[jj])
          else:
            fp.write('%1.15e\t' % p_jj[jj])

      toc = timeit.default_timer()
      print "timer = %1.14e"%(toc-tic)
      fp.close()

    if(return_picks):
      # TODO: hack: better to never create file, or have file generation be a separte routine
      picks = np.loadtxt(filename)
      os.remove(filename)
      return picks

    else:
      raise Exception("sampling type unknown")


def remap_parameters(file_in,file_out,param_dict_map={}):
  """ param_dict_map -- param_dict_map[i] is a function from 
                        an n-dimensional vector of parameter values
                        to the new ith parameter.
      file_in --  file where each row is P=(p[0],p[1],...,p[n-1])
      file_out -- each row is now 
                  (param_dict_map["0"](P),...,param_dict_map["n-1"](P)"""


  old_params = np.loadtxt(file_in)

  total_params  = old_params.shape[0]
  parameter_dim = old_params.shape[1]

  print "parameter dimensionality = ",parameter_dim
  print "total parameters = ",total_params

  fp = open(file_out,'w')

  for ii in range(total_params):

    pii = old_params[ii,:]

    for jj in range(parameter_dim):

      # apply parameter mapping if requested
      if param_dict_map.has_key(str(jj)):
        p_jj = param_dict_map[str(jj)](pii)
      else:
        p_jj = pii[jj]

      if(jj == parameter_dim-1):
        fp.write('%1.15e\n' % p_jj)
      else:
        fp.write('%1.15e\t' % p_jj)

  fp.close()


if __name__=="__main__":

    ### typical usage: 3 random sample files are desired...
    ### enter base name as first argument and 3 as second

    filename = str(sys.argv[1])
    files    = int(sys.argv[2])
    for ii in range(files):
        generate_sampling(filename+str(ii))


