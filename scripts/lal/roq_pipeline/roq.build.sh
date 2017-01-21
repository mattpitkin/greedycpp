#PBS -l nodes=1
#PBS -l walltime=2:30:00
#PBS -o seobv2_roq5_quad.out
#PBS -e seobv2_roq5_quad.err
#PBS -A sxs
#PBS -d .
#PBS -S /bin/bash
#PBS -N build.basis

let "processes=1*12"
bindir="/home/sfield/greedycpp/bin" # greedy executable here

### USER SETTINGS HERE ###

# phenomPv2 runs
#cfgpath="/panfs/ds08/sxs/sfield/greedy_phenomPv2/phenomp_PlusCross_parts_4s/greedyTS.cfg"
#cfgpath="/panfs/ds08/sxs/sfield/greedy_phenomPv2/phenomp_Quad_parts_128s/run_settings.cfg"

# Lackey tidal and SEOB runs
source /home/sfield/zwicky-env/SetupGreedyICC.sh
module add gsl/1.16
cfgpath="default_settings.cfg"
#cfgpath="/panfs/ds08/sxs/sfield/greedy_seob/production_runs/roq_5/df128/run_settings.cfg"
#cfgpath="/panfs/ds08/sxs/sfield/greedy_tidalEOB/lea_review/run_settings.cfg"
#cfgpath="/panfs/ds08/sxs/sfield/greedy_seob/production_runs/roq_5/df128/ts_quad/run_settings.cfg"

##########################

### DERIVED SETTINGS HERE ###
outdir=`grep "output_dir" ${cfgpath} | grep -oP '(?<=").+?(?=")'` # folder where data will be place

#module add openmpi/1.4.1-intel
#module add openmpi/1.4.3-gcc-static
echo `env`

echo "processes = " ${processes}
echo "outdir is " ${outdir}


NPROCS=`wc -l $PBS_NODEFILE |gawk '//{print $1}'`
echo "* running on $NPROCS cores..."
echo Master process running on `hostname`
echo Directory is `pwd`
echo Starting execution at `date`

# Display libraries; this will be written to your stdout file.
/usr/bin/ldd $greedyexe
cat $PBS_NODEFILE


# NOTE: comment out APP with three "###"
### basic mpi, 1 process per core ###
APP="mpirun -np +${processes} ${bindir}/greedy_mpi ${cfgpath}"
#APP="${bindir}/greedy ${cfgpath}"
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

#wait


