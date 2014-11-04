''' Generate parameter sampling, write these values to file. The resulting
    text file can be read as input to greedy.cpp through proper specification
    in the configuration file.'''

import numpy as np
import math
import sys

### output file name ###

def generate_sampling(filename):

    ### sampling strategies (linear, log10, ln) ###
    sample_type    = "rand" # "rand" or "deterministc"

    ### setup parameter intervals ###
    params_low  = [1.0,1.0]  # lower interval of each parameter
    params_high = [3.0,3.0]  # params_num[i] is the number of samplings in the interval [param_low[i],param_high[i]]
 
    ### setup for deterministic sampling ###
    param_sampling = "ln"
    scaling_factor = 20.0 # default should be 1. higher values densely sample lower range of parameter interval
    params_num  = [50,50] # deterministic: upper interval of each parameter

    ### setup for random sampling ###
    total_picks = 2000     # random: this many draws from interval

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

      ### output MyTS.txt here -- tensor product parameter grid ###
      fp = open(filename,'w')
      for ii in range(np.size(p1)):
        for jj in range(np.size(p2)):
          fp.write('%1.15e\t' % p1[ii])
          fp.write('%1.15e\n' % p2[jj])

      fp.close()

    elif( sample_type is "rand"):

      fp = open(filename,'w')

      for ii in range(total_picks):
        for jj in range(len(params_high)):

          p_jj = (params_high[jj] - params_low[jj]) * np.random.random_sample((1,)) + params_low[jj]

          if(jj == len(params_high)-1):
            fp.write('%1.15e\n' % p_jj)
          else:
            fp.write('%1.15e\t' % p_jj)

    else:
      raise Exception("sampling type unknown")



if __name__=="__main__":

    filename = str(sys.argv[1])
    generate_sampling(filename)
