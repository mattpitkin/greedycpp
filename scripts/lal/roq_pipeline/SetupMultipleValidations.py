#!/usr/bin/env python
import os
from string import Template

# if multiple roq "bands" (sets VALINTERVALS, rundir, and random points directory)
# BASISDIR still needs to be set by hand
roq_indx='5'

# name of directory for validation run setup 
#rundir_name='rundir'
rundir_name='rundir'+roq_indx

Simulations = 320
#Simulations = 320

replacement_dict = {}
replacement_dict['NODES']='1'

#replacement_dict['WALLTIME']='10:30:00' # 40 - 4096, 400000, ++, xx, +x
#replacement_dict['WALLTIME']='00:30:00' # 40 - 4096, 400000, ++, xx, +x

# NOTE: make sure to use > 2 hours (otherwise the debug queue will get clogged up)
replacement_dict['WALLTIME']='10:00:00' # 40 - 4096, 400000, +, x

# WARNING - NOTE - READ: 
#  location of basis and eim files - MUST END WITH "/" #
#replacement_dict['BASISDIR']='/panfs/ds08/sxs/sfield/greedy_phenomPv2/phenomp_Quad_parts_128s/'
#replacement_dict['BASISDIR']='/panfs/ds08/sxs/sfield/greedy_tidalEOB/Lackey_qaud_20hz/'
#replacement_dict['BASISDIR']='/panfs/ds08/sxs/sfield/greedy_tidalEOB/Lackey_qaud_20hz/df128/'
#replacement_dict['BASISDIR']='/panfs/ds08/sxs/sfield/greedy_tidalEOB/Lackey_PlusCross_20hz/df128/'
#replacement_dict['BASISDIR']='/panfs/ds08/sxs/sfield/greedy_seob/production_runs/roq_3/df32/ts_quad/'
replacement_dict['BASISDIR']='/panfs/ds08/sxs/sfield/greedy_seob/production_runs/roq_5/df128/ts_quad/'

# file containing intervals on which to validate -- can be left blank? #
#replacement_dict['VALINTERVALS']="/home/sfield/misc_notebooks_codes/phenomP_cfg_scripts/data/phenomPv2/non_kerr_violating_run/128s/validation_range.dat"
#replacement_dict['VALINTERVALS']="junk.junk"
replacement_dict['VALINTERVALS']="/panfs/ds08/sxs/sfield/greedy_seob/production_runs/roq_"+roq_indx+"/mass_partition_"+roq_indx+".input"

# where validation results will go - existing files of the same name overwritten #
replacement_dict['VALRESULTSDIR']="validations_full_range"

# validation cfg file name #
replacement_dict['VALCFG']="validations_setup.cfg"

# script needs to take numberOfSamples, outputFilename, validation intervals (given above)
#replacement_dict['SAMPLER_SCRIPT']="/home/sfield/misc_notebooks_codes/phenomP_cfg_scripts/helper_scripts/sample_phenomP.py"
#replacement_dict['SAMPLER_SCRIPT']="/home/sfield/tidalroq/src/scripts/makeTS_power_boundary_component_masses.py"
replacement_dict['SAMPLER_SCRIPT']="/home/sfield/misc_notebooks_codes/greedycpp_helpers/examples/seob_roq/generate_seob_points.py"

# location of validation executable #
replacement_dict['BINDIR']='/home/sfield/greedycpp/bin'

### randfile: file of random parameter values -- create a new one or point to existing one
replacement_dict['RANDFILE']="/panfs/ds08/sxs/sfield/tmp"+roq_indx+"/test.txt"
#replacement_dict['RANDFILE']="/panfs/ds08/sxs/sfield/tmp/test.txt"

replacement_dict['MODE_MOD_SCRIPT']="/home/sfield/misc_notebooks_codes/greedycpp_helpers/AddReplaceModeDimension.py"

#replacement_dict['SAMPLES']='1000'
#replacement_dict['SAMPLES']='25000'
replacement_dict['SAMPLES']='30000' # times 3 (about 1 in 4 picks are retained) for B20 validation
#replacement_dict['MODES']="+ x"
#replacement_dict['MODES']="+" # spin-aligned systems have hx = (complex number) hx
#replacement_dict['MODES']="++ xx +x"
replacement_dict['MODES']="++ xx h2" # anntena patterns (1,0), (0,1), (1,1)

#######################################################################
val_template = """
#PBS -l nodes=$NODES
#PBS -l walltime=$WALLTIME
#PBS -A sxs
#PBS -o $BASISDIR$VALRESULTSDIR$SIMNUM/val.out.$SIMNUM
#PBS -e $BASISDIR$VALRESULTSDIR$SIMNUM/val.err.$SIMNUM
#PBS -d .
#PBS -S /bin/bash
#PBS -N validation.tidalEOB
#PBS -m ea
#PBS -M sef74@cornell.edu

# Lackey tidal EOB runs need this
source /home/sfield/zwicky-env/SetupGreedyICC.sh
module add gsl/1.16

## run the sampler ##
# phenomPv2 format
#${SAMPLER_SCRIPT} --p $SAMPLES --f $RANDFILE$SIMNUM --v $VALINTERVALS
# Lackey Tidal EOB format
#${SAMPLER_SCRIPT} --N $SAMPLES --f $RANDFILE$SIMNUM --sampler random --full
# seob format
${SAMPLER_SCRIPT} --N $SAMPLES --o $RANDFILE$SIMNUM --input $VALINTERVALS --sampler random 

## add modes if model is "all parts" ##
${MODE_MOD_SCRIPT} --f $RANDFILE$SIMNUM --m $MODES

NPROCS=`wc -l $PBS_NODEFILE |gawk '//{print $1}'`
echo "* running on $NPROCS cores..."
echo Master process running on `hostname`
echo Directory is `pwd`
echo Starting execution at `date`

cd $BINDIR

# NOTE: comment out APP with three "###"
APP="./verifyOMP $BASISDIR$VALCFG $BASISDIR gsl $RANDFILE$SIMNUM $VALRESULTSDIR$SIMNUM"
echo " "
echo $APP
echo " "
$APP


echo Ending execution at `date`
echo "* done"
"""

#######################################################################

# setup the run directory #
rwd = os.getcwd()+'/'+rundir_name
os.makedirs(rwd)

jobs = open(rwd+'/submit_all_jobs.sh','w')
merge_script = open(rwd+'/merge_validations.sh','w')
DESTDIR = replacement_dict['BASISDIR']+replacement_dict['VALRESULTSDIR']
merge_script.write('mkdir '+DESTDIR+'\n')

for i in range(Simulations):

  replacement_dict['SIMNUM'] = str(i)

  run_dir = rwd+"/run_dir"+str(i)
  os.makedirs(run_dir)
  template = Template(val_template)
  submission_string = template.safe_substitute(replacement_dict)

  fp = file(run_dir+'/roq.validations.sh','w')
  fp.write(submission_string)
  fp.close()

  # add job to job submission script
  jobs.write('qsub '+run_dir+'/roq.validations.sh\n')

  # add validation merger directions to script
  merge_script.write('mv '+DESTDIR+str(i)+'/test.txt'+str(i)+' '+DESTDIR+'\n') # TODO: hardcoded RANDFILE
  merge_script.write('mv '+DESTDIR+str(i)+'/bad_points_test.txt'+str(i)+' '+DESTDIR+'\n') # TODO: hardcoded RANDFILE
  merge_script.write('mv '+DESTDIR+str(i)+'/val.out.'+str(i)+' '+DESTDIR+'\n')
  merge_script.write('mv '+DESTDIR+str(i)+'/val.err.'+str(i)+' '+DESTDIR+'\n')
  if i==0:
    merge_script.write('mv '+DESTDIR+str(i)+'/validation_run.cfg'+' '+DESTDIR+'\n')

jobs.close()

# collect list of bad points (points with large errors).
merge_script.write('cat '+DESTDIR+'/bad_points_test.txt* > '+DESTDIR+'/ts_bad_points.txt\n')
merge_script.write('rm '+DESTDIR+'/bad_points_test.txt*\n')

# Combined greedy points + bad points = next training set (iteration = it)
next_it_dir = replacement_dict['BASISDIR']+'/it/'
merge_script.write('mkdir '+next_it_dir+'\n')
merge_script.write('cat '+replacement_dict['BASISDIR']+'/GreedyPoints.txt '+DESTDIR+'/ts_bad_points.txt > '+next_it_dir+'TS_IT.txt\n')
merge_script.close()

os.chmod(rwd+'/submit_all_jobs.sh',0700) # make job submission script executable
os.chmod(rwd+'/merge_validations.sh',0700) # make job submission script executable

