# output files for greedy runs #
#outdir="/panfs/ds06/sxs/sfield/greedy_phenomP/phenomp_plus_20_40/" # location of all output files
#outdir="/panfs/ds06/sxs/sfield/greedy_phenomP/phenomp_cross_20_40/" # location of all output files
#outdir="/panfs/ds06/sxs/sfield/greedy_phenomP/phenomp_plus_15_20/"
#outdir="/panfs/ds06/sxs/sfield/greedy_phenomP/phenomp_cross_15_20/"
#outdir="/panfs/ds06/sxs/sfield/greedy_phenomP/phenomp_PlusCross_parts_40_4096/"
#outdir="/panfs/ds06/sxs/sfield/greedy_phenomP/phenomp_PlusCross_parts_15_4096/"
outdir="/panfs/ds06/sxs/sfield/greedy_phenomP/phenomp_PlusCross_parts_40_4096_iteration1/"


# random sample output directory (if multiple runs) #
#randdir="plus_validation_full_range"
#randdir="plus_validation"
randdir="validations_both_parts_full_range"

# validation cfg file name #
#valcfg="run_settings.cfg"
valcfg="validations_setup.cfg"
#valcfg="validations_plus.cfg"
#valcfg="validations_cross.cfg"


run_dir='/home/sfield/greedycpp/bin'

PathToValidations=${outdir}${randdir}

echo ${PathToValidations}

#PBS -l nodes=1
#PBS -l walltime=00:30:00
#PBS -A sxs
#PBS -o ${PathToValidations}/val.out
#PBS -e ${PathToValidations}/val.err
#PBS -d .
#PBS -S /bin/bash
#PBS -N validation.phenomP

#threads=12 #12(1) 108(9) 144(12) 192(16) 264(22) 360(30) 480(40) 600(50) 

### randfile -- create a new one
randfile="/panfs/ds06/sxs/sfield/tmp/test.txt"

# random sample files -- existing files #
#randfile="Rand_phenom_GP_cross_20_40_1e6_smalltest.txt" ## test file
#randfile="Rand_phenom_GP_cross_20_40_1e6_set1.txt"
#randfile="Rand_phenom_GP_cross_or_plus_10_20_1e6.txt"
#randfile="Rand_phenom_GP_full_parameter_range_1e6.txt"


cat $PE_HOSTFILE > hostfile.txt
awk '{$1=$1; print $1}' hostfile.txt > machfile


## run the sampler ##
/home/sfield/misc_notebooks_codes/phenomP_cfg_scripts/helper_scripts/sample_phenomP.py --p 1100 --f ${randfile}

## add modes if model is "all parts" ##
/home/sfield/misc_notebooks_codes/phenomP_cfg_scripts/helper_scripts/AddModeDimension.py --f ${randfile} --m + x


NPROCS=`wc -l $PBS_NODEFILE |gawk '//{print $1}'`
echo "* running on $NPROCS cores..."
echo Master process running on `hostname`
echo Directory is `pwd`
echo Starting execution at `date`

cd $run_dir

# NOTE: comment out APP with three "###"
APP="mpirun -np 1 ./verifyOMP ${outdir}${valcfg} ${outdir} gsl ${randfile} ${randdir}"
echo " "
echo $APP
echo " "
$APP


echo Ending execution at `date`
echo "* done"
