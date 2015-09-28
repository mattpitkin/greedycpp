

### USER SETTINGS HERE ###
cfgpath="/panfs/ds06/sxs/sfield/greedy_phenomP/phenomp_PlusCross_parts_40_4096_iteration1/run_settings.cfg"
bindir="/home/sfield/greedycpp/bin" # greedy executable here

#PBS -l nodes=4
let "threads=4*12"

##########################

### DERIVED SETTINGS HERE ###
outdir=`grep "output_dir" ${cfgpath} | grep -oP '(?<=").+?(?=")'` # folder where data will be place


#PBS -l walltime=03:00:00
#PBS -A sxs
#PBS -o ${outdir}/greedybuild.out
#PBS -e ${outdir}/greedybuild.err
#PBS -d .
#PBS -S /bin/bash
#PBS -N build.basis

echo "threads = " ${threads}
echo "outdir is " ${outdir}


NPROCS=`wc -l $PBS_NODEFILE |gawk '//{print $1}'`
echo "* running on $NPROCS cores..."
echo Master process running on `hostname`
echo Directory is `pwd`
echo Starting execution at `date`

#cd $run_dir

# NOTE: comment out APP with three "###"
### basic mpi, 1 process per core ###
APP="mpirun -np +${threads} ${bindir}/greedy_mpi ${cfgpath}"
echo " "
echo $APP
echo " "
$APP

echo Ending basis building at `date`

# NOTE: second input to eim, the basis+quadrature directory, must end in '/'
APP="${bindir}/eim ${cfgpath} ${outdir}/ gsl"
echo " "
echo $APP
echo " "
$APP

echo Ending execution at `date`
echo "* done"
