''' Generate 2D parameter sampling, write to file MyTS.txt. The resulting
    text file can be read as input to greedy.cpp through proper specification
    in the configuration file.'''

import numpy as np
import math

### sampling strategies (linear, log10, ln) ###
param_sampling = "ln";
scaling_factor = 20.0 # default should be 1. higher values densely sample lower range of parameter interval

### setup parameter intervals ###
params_num  = [500,500]  # params_num[i] is the number of samplings in the interval [param_low[i],param_high[i]]
params_low  = [1.0,1.0]  # lower interval of each parameter
params_high = [4.0,4.0]  # upper interval of each parameter


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
fp = open('MyTS.txt','w')
for ii in range(np.size(p1)):
  for jj in range(np.size(p2)):
    fp.write('%1.15e\t' % p1[ii])
    fp.write('%1.15e\n' % p2[jj])

fp.close()
