''' Generate 2D parameter sampling, write to file MyTS.txt. The resulting
    text file can be read as input to greedy.cpp through proper specification
    in the configuration file.'''

import numpy as np

### setup parameter intervals ###
params_num  = [30,30]    # params_num[i] is the number of samplings in the interval [param_low[i],param_high[i]]
params_low  = [1.0,1.0]  # lower interval of each parameter
params_high = [3.0,3.0]  # upper interval of each parameter

### sample the parameter space ###
# linear spacing of p1 and p2
p1 = np.linspace(params_low[0],params_high[0],params_num[0])
p2 = np.linspace(params_low[1],params_high[1],params_num[1])
# logarithmic spacing, clusters more points near lower end of the interval
#p1 = np.power(params_low[0]*(params_high[0]/params_low[0]),np.linspace(0,1,params_num[0]))
#p2 = np.power(params_low[1]*(params_high[1]/params_low[1]),np.linspace(0,1,params_num[1]))

### output MyTS.txt here -- tensor product parameter grid ###
fp = open('MyTS.txt','w')
for ii in range(np.size(p1)):
  for jj in range(np.size(p2)):
    fp.write('%1.15e\t' % p1[ii])
    fp.write('%1.15e\n' % p2[jj])

fp.close()
