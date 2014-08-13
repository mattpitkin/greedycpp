import numpy as np

### setup parameter intervals ###
params_num  = [30,30]    # params_num[i] is the number of samplings in the interval [param_low[i],param_high[i]]
params_low  = [1.0,1.0]  # lower interval of each parameter (m1 and m2)
params_high = [3.0,3.0]  # upper interval of each parameter (m1 and m2)

### sample the parameter space ###
p1 = np.linspace(params_low[0],params_high[0],params_num[0])
p2 = np.linspace(params_low[1],params_high[1],params_num[1])
#p1 = np.power(params_low[0]*(params_high[0]/params_low[0]),np.linspace(0,1,params_num[0]))
#p2 = np.power(params_low[1]*(params_high[1]/params_low[1]),np.linspace(0,1,params_num[1]))

### output MyTS.txt here ###
fp = open('MyTS.txt','w')
for ii in range(np.size(p1)):
	for jj in range(np.size(p2)):
		fp.write('%1.15e\t' % p1[ii])
		fp.write('%1.15e\n' % p2[jj])

fp.close()
