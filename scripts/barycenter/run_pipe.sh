#!/bin/sh

# script to run greedycpp, then get the empirical interpolator, then run the validator, and finally enrich the original training set!
# (see e.g. https://bitbucket.org/sfield83/greedycpp/wiki/features/Generic.md)

source /home/matthew/repositories/greedycpp/script/barycenter/setup.sh

# set the executables
genscript=/home/matthew/repositories/greedycpp/scripts/barycenter/generate_trainingpoints.py # script to generate set of RA and dec sky positions as training points
greedy=/home/matthew/repositories/greedycpp/bin/greedyOMP  # greedycpp executable
eim=/home/matthew/repositories/greedycpp/bin/eim           # empirical interpolator exectuable
verify=/home/matthew/repositories/greedycpp/bin/verifyBarycenter  # verification executable
#verify=/home/matthew/repositories/greedycpp/bin/verifyOMP  # verification executable

rundir=/home/matthew/testing/redordbar/
if [ ! -d "$rundir" ]; then
  mkdir -p $rundir
fi
cd $rundir

# training set output file
tspoints=$rundir/TS_points.txt

# generate a training set of sky positions
nsky=10000
$genscript -N $nsky -o $tspoints

# set the greedycpp configuration file
configfilename=test.cfg
configfile=$rundir/$configfilename

# run greedy cpp
$greedy $configfile

# run the EIM code
$eim $configfile $rundir gsl

# set the validation directory
valdirname=validation
valdir=${rundir}/${valdirname}
mkdir -p $valdir
valpoints=${valdir}/validation_points.txt

# create the validation points
$genscript -N $nsky -o $valpoints

$verify ${rundir}/validations_setup.cfg $rundir gsl $valpoints $valdirname

# run enrichment with a concatenation of the originally output greedy points and the "bad" points
enrichdir1=${rundir}/enrich1/
mkdir -p $enrichdir1
cat $rundir/GreedyPoints.txt $valdir/bad_points_validation_points.txt > $enrichdir1/TS_points.txt
cp $configfile $enrichdir1
configfilenew=$enrichdir1/$configfilename

cd $enrichdir1

# re-run greedy algorithm
$greedy $configfilenew

# re-run the EIM code
$eim $configfile $enrichdir1 gsl

valdirnew=${enrichdir1}/${valdirname}
mkdir -p $valdirnew
valpointsnew=${valdirnew}/validation_points.txt

# create the validation points
$genscript -N $nsky -o $valpointsnew

$verify ${enrichdir1}/validations_setup.cfg $enrichdir1 gsl $valpointsnew $valdirname

