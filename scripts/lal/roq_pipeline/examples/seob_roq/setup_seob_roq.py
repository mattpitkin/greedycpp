#!/usr/bin/env python
''' Setup a handful of ROQ traning intervals using *.input files in PARTION_DIR '''

import subprocess, os, json, fnmatch
from BatchSubmission.SimulationDataTools import fancy_glob
import numpy as np



# USER SETTINGS HERE -- common to all ROQ training intervals #
roq_sim = {
#    "ts sampling": "boundary",
    "ts sampling": "full",
    "MODEL": "SEOBNRv2_ROM_DoubleSpin_HI_all_parts",
    "modes": "+ x",
    "sampler script": "/home/sfield/misc_notebooks_codes/greedycpp_helpers/examples/seob_roq/generate_seob_points.py",
    "partitions": []
}

partition_dir = "/panfs/ds08/sxs/sfield/greedy_seob/production_runs/partitions_roq/"
#partition_dir = "/panfs/ds08/sxs/sfield/greedy_seob/test_runs/PipelineROQ/partitions_roq/"

partition_files = fancy_glob(partition_dir+'mass_partition_*.input',sort=True)
print len(partition_files)

# common frequency interval #
#flow = 20
#fup  = 4096
# or frequency interval for each mass partition 
try:
  f_intervals = np.loadtxt(partition_dir+'frequency_intervals.input')
except IOError:
  f_intervals = np.zeros([len(partition_files),2])
  f_intervals[:,0] = flow
  f_intervals[:,1] = fup

print f_intervals



################################################################

for i, f1 in enumerate(partition_files):

  f = f1[1]
  flow = f_intervals[i,0]
  fup  = f_intervals[i,1]

  print("File = %s with flow = %f, fup = %f"%(f,flow,fup))

  sample_intervals = partition_dir+f

  subprocess.call('/home/sfield/misc_notebooks_codes/greedycpp_helpers/examples/seob_roq/generate_seob_quadrature.py -l '+\
                  str(flow)+' -u '+str(fup)+' -i '+f, shell=True)
  #subprocess.call('/home/sfield/misc_notebooks_codes/greedycpp_helpers/examples/seob_roq/generate_seob_quadrature.py -l '+\
  #                str(flow)+' -u '+str(fup)+' -i '+f+' --mlow 1.2', shell=True)

  roq_sim['partitions'].append({
            "quad": "adaptive",
            "quad_nodes_file": f.split('.')[0]+'_Nodes.dat',
            "quad_weight_file": f.split('.')[0]+'_Weights.dat',
            "sample intervals": f
        })


with open('roq_info.json', 'w') as fp:
  json.dump(roq_sim, fp, sort_keys=True,indent=4, separators=(',', ': '))
