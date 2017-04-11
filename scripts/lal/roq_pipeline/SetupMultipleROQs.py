#!/usr/bin/env python
import subprocess, os, json, shutil
from string import Template


### Read in setup from json file ###
json_file = 'roq_info.json'
with open(json_file) as fp:
  roq_simulations = json.load(fp)
print roq_simulations


ROQ_SCRIPTS='/home/sfield/greedycpp/scripts/lal/roq_pipeline/'

#MODEL = "Lackey"

#SAMPLER_SCRIPT = "/home/sfield/misc_notebooks_codes/phenomP_cfg_scripts/helper_scripts/sample_phenomP.py"
#SAMPLE_INTERVALS = "/home/sfield/misc_notebooks_codes/phenomP_cfg_scripts/data/phenomPv2/non_kerr_violating_run/128s/validation_range.dat"

### setup multiple runs from dictionary (either specify below or read from JSON file) ###
replacement_dict = {}
#"LackeyTidal2013_SEOBNRv2_ROM_HI_all_parts";
#replacement_dict['partitions'] = [{'quad':'adaptive',
#'quad_nodes_file':'\"/panfs/ds08/sxs/sfield/greedy_tidalEOB/input/QuadratureRule_20Hz_smooth_fine.dat\";',
#'quad_weight_file':'\"/panfs/ds08/sxs/sfield/greedy_tidalEOB/input/MyWeights_20Hz_smooth_fine.dat\";'
#},
#{'quad':'uniform',
#'x_min':'20.0',
#'x_max':'4096.0',
#'quad_points':'521729'
#},
#]


#roq_simulations = {
#'MODEL':"SEOBNRv2_ROM_DoubleSpin_HI_all_parts",
#'modes': '+',
#'sampler script': ROQ_SCRIPTS+"generate_seob_points.py",
#'partitions':[
#{'quad':'adaptive',
#'sample intervals': ROQ_SCRIPTS+"seob_sample.input",
#'quad_nodes_file':'/home/mpuer/projects/greedycpp/runs_MP/QuadratureRule_SEOBNRv2_20Hz_smooth_fine.dat',
#'quad_weight_file':'/home/mpuer/projects/greedycpp/runs_MP/MyWeights_SEOBNRv2_20Hz_smooth_fine.dat',
#}]
#}

#with open('roq_info.json', 'w') as fp:
#  json.dump(roq_simulations, fp, sort_keys=True,indent=4, separators=(',', ': '))

### Nothing to change below this line ###

MODE_MOD_SCRIPT = ROQ_SCRIPTS+"AddReplaceModeDimension.py"

#######################################################################

### Model specific parts of input file ###
cfg_seob_rom = """
// parameter domain and its sampling 
param_dim      = 5;                       // (int) number of paramteric dimensions (currently supports 2)
load_from_file = true;                   // (bool) load training points from file
// params = {m1,m2,chi1,chi2}
p1_scale       = 1.0; // BH1 mass in Solar masses
p2_scale       = 1.0; // BH2 mass in Solar masses 
p3_scale       = 1.0; // dimensionless spin of BH1
p4_scale       = 1.0; // dimensionless spin of BH2
p5_scale       = 1.0; // GW's polarization type for all_parts model

model_name = "SEOBNRv2_ROM_DoubleSpin_HI_all_parts";
"""

cfg_lackey = """
// parameter domain and its sampling 
param_dim      = 5;                       // (int) number of paramteric dimensions
load_from_file = true;                   // (bool) load training points from file
// params = {mBH,mNS,chiBH,Lambda}
p1_scale       = 1.0; // BH mass in Solar masses
p2_scale       = 1.0; // NS mass in Solar masses 
p3_scale       = 1.0; // dimensionless spin of BH
p4_scale       = 1.0; // tidal Lambda
p5_scale       = 1.0; // GW's polarization type for all_parts model

model_name = "LackeyTidal2013_SEOBNRv2_ROM_HI_all_parts";
"""

### generic portions of input file ###
cfg_general = """
// greedy algorithm information
seed           = 0;                     // (int) greedy algorithm global index seed
tol            = 1e-11;                 // (double) greedy algorithm tolerance achieving \| app \|^2
max_RB         = 3000;                 // (int) estimated number of RB (reasonable upper bound)

output_dir         = "${full_path}";
output_data_format = "bin";

ts_file             = "${full_path}/TS.txt"

export_tspace = false;
"""

quadrature_template_adaptive = """// inner product information
quad_type       = 2;                  // (int) 0 = LGL, 1 = Reimman sum, 2 = user-defined (via input files)
quad_nodes_file = "${quad_nodes_file}";
num_weight_file = "${quad_weight_file}";
weighted_inner  = false;              // (bool) whether to include a weight W(x): \int W(x) f(x) g(x)
"""

quadrature_template_uniform = """
quad_type = 1;
x_min = ${x_min};
x_max = ${x_max};
quad_points = ${quad_points};
weighted_inner  = false;
"""

#######################################################################

# setup the run directory #
rwd = os.getcwd()+'/'

MODES = roq_simulations['modes']
MODEL = roq_simulations['MODEL']
SAMPLER_SCRIPT = roq_simulations['sampler script']
SAMPLING = roq_simulations["ts sampling"]

shutil.copy(ROQ_SCRIPTS+'SetupMultipleValidations.py',rwd)

for i, partition in enumerate(roq_simulations['partitions']):

  print partition

  rundir                 = 'roq_'+str(i)
  full_path              = rwd+rundir
  partition['full_path'] = full_path
  SAMPLE_INTERVALS = partition['sample intervals']


  if partition['quad']=='adaptive':
    cfg = quadrature_template_adaptive + cfg_general
  elif  partition['quad']=='uniform':
    cfg = quadrature_template_uniform + cfg_general
  else:
    raise ValueError

  if MODEL == "LackeyTidal2013_SEOBNRv2_ROM_HI_all_parts":
    cfg += cfg_lackey
  elif MODEL == "SEOBNRv2_ROM_DoubleSpin_HI_all_parts":
    cfg += cfg_seob_rom
  else:
    raise ValueError

  os.makedirs(full_path)
  template = Template(cfg)
  submission_string = template.safe_substitute(partition)

  # write cfg file
  fp = file(full_path+'/default_settings.cfg','w')
  fp.write(submission_string)
  fp.close()

  # copy files into local run directories
  shutil.copy(SAMPLE_INTERVALS,full_path)
  shutil.copy(partition['quad_nodes_file'],full_path)
  shutil.copy(partition['quad_weight_file'],full_path)
  shutil.copy(ROQ_SCRIPTS+'roq.build.sh',full_path)

  if SAMPLING == "boundary":
    subprocess.call(SAMPLER_SCRIPT+' --no-full -o'+rundir+'/TS.txt --input '+SAMPLE_INTERVALS, shell=True)
  elif SAMPLING == "full":
    subprocess.call(SAMPLER_SCRIPT+' --f -o '+rundir+'/TS.txt --input '+SAMPLE_INTERVALS, shell=True)
  else:
    raise valueError

  subprocess.call(MODE_MOD_SCRIPT+' --f '+rundir+'/TS.txt --m '+MODES, shell=True)


# clean up aux file
os.remove('junk.txt')
