#!/bin/bash
#
# queue options:
#    runtime: -l h_rt=hr:min:sec
#$ -l h_rt=24:00:00
#    number of procs (12 per node): -pe [type] [num]
#$ -pe orte 24
#$ -N my_job
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o test.out

# this environment variable is set as the directory where the job was submitted
#cd $SGE_O_WORKDIR

echo "* running on $NSLOTS cores..."
echo Master process running on `hostname`
echo Directory is `pwd`
echo Starting execution at `date`

### to generate a quadrature rule ###
#python generate_quadrature_rule.py


# SEOB #
#rundir="/home/sfield/greedy-seob/models" #cd into this directory
#bindir="/home/sfield/greedy-seob" # greedy executable here
#cfgpath="/home/sfield/greedy-seob/SEOB/EOB_SS.cfg"


# greedycpp #
rundir="/home/sfield/greedycpp" #cd into this directory
bindir="/home/sfield/greedycpp" # greedy executable here
outdir="/data/sfield/my_output_GW40hzTEST" # see *.cfg file --> machfile and output file copied here
#cfgpath="/home/sfield/greedycpp/examples/test1.cfg"
cfgpath="/home/sfield/greedycpp/examples/TaylorF2_BNS_40Hz/test1.cfg"

### these will get copied into outdir ###
machpath="/home/sfield/machfile"
outfile="/home/sfield/test.out"


cd $rundir


### basic mpi, 1 process per core ###
mpirun -np 132 -report-bindings $bindir/greedy_mpi $cfgpath



### large orthogonalization run: 1 process for ortho (rank=0) node, 1 proces/core for workers on other nodes  ###

# here we (i) grab the hostfile, put 1 node per line, enforce 12 processes per node and 1 process on the head node 
#cat $PE_HOSTFILE > hostfile.txt
#awk '{$1=$1; print $1}' hostfile.txt > machfile
#sed -i "s/$/ slots=12 max-slots=12/" machfile
#sed -i "1s/ slots=12 max-slots=12/ slots=1 max-slots=1/" machfile
#let "procs = $NSLOTS - 11" # orthogonalization node (rank 0) has 1 process only 

# hostfile type 1
#mpirun -np $procs --hostfile $machpath -report-bindings $bindir/greedy_mpi $cfgpath

###1 node for ortho, 1 node for worker process and openMP threads for ts sweep ###
# Set OMP_NUM_THREADS to be less than cores per node (e.g. 11)
# hostfile type 2
#mpirun -np $procs --hostfile $machpath -report-bindings $bindir/greedy_mpi $cfgpath



### verification jobs ###
#./../verify /data/sfield/SEOB_MM99_iteration2_output/EOB_SS.cfg /data/sfield/SEOB_MM99_iteration2_output/ /home/sfield/greedy-seob/scripts/rand.txt13
#./verify /data/sfield/my_output_GW40hz/test1.cfg /data/sfield/my_output_GW40hz/ /home/sfield/greedycpp/scripts/random8.txt 
#./verify my_output_test1/test1.cfg my_output_test1/ scripts/random.txt


echo Ending execution at `date`
echo "* done"

cp test.out $outdir
cp $machpath $outdir


