#!/usr/bin/env python

"""
Create a pipeline, run as a Condor DAG, to perform the ROM creation, EIM calculation, validation and (if required) enrichment of the ROM,
for barycentring time delays.
"""

from __future__ import division, print_function

import os
import sys
import argparse

"""
things to add:
 - subfile for running greedyOMP that includes "request_memory = X" and "machine_count = Y" for running on multiple cores
 - create the .cfg file for running the code
 - check whether to validate
 - check whether to enrich (and set a maximum number of enrichment steps)
"""

configtemplate = """// Configuration file for running greedycpp with the 'Barycenter' model

// start time, end time and number of time points
x_min = {starttime};     // GPS start time
x_max = {endtime};       // GPS end time
quad_points = {npoints}; // number of time points

// model information
model_name = "Barycenter_{detector}_{ephemeris}_{timeunits}";

// inner product information
quad_type = 1;
weighted_inner = false;

param_dim = 2;
load_from_file = true; // load training points (RA and dec values) from file
p1_scale = 1.0;        // scaling for right ascension
p2_scale = 1.0;        // scaling for declination

ts_file = "{tsfile}";  // location of training points file

// greedy algoritm information
seed = 0;
tol = {tolerance};  // greedy algorithm tolerance
max_RB = {maxrb};   // estimated upper bound in number of reduced bases

// output information
outputdir = "{outdir}";      // output directory
output_data_format = "both"; // output both in binary and text format
"""

parser = argparse.ArgumentParser( )
parser.add_argument("--rundir", dest="rundir", required=True, help="Set the run directory for the analysis")
parser.add_argument("--execdir", dest="execdir", required=True, help="Set the directory containing the greedycpp executables")
parser.add_argument("--ntraining", dest="ntraining", type=int, default=5000, help="The number of training points for the greedy algorithm [default: %(default)s]")
parser.add_argument("--starttime", dest="starttime", type=float, default=900000000.0, help="GPS start time of templates [default: %(default)s]")
parser.add_argument("--endtime", dest="endtime", type=float, default=900863400.0, help="GPS end time of templates [default: %(default)s]")
parser.add_argument("--npoints", dest="npoints", type=int, default=1440, help="The number of time points for templates [default: %(default)s]")
parser.add_argument("--detector", dest="detector", default="H1", help="Gravitational wave detector to use [default: %(default)s]")
parser.add_argument("--ephemeris", dest="ephem", default="DE405", help="JPL solar system ephemeris to use (DE200, DE405, DE414, or DE421) [default: %(default)s]")
parser.add_argument("--units", dest="units", default="TCB", help="Time units to use (TCB or TDB) [default: %(default)s]")
parser.add_argument("--tolerance", dest="tol", type=float, default=5e-12, help="The tolerance for producing the reduced basis [default: %(default)s]")
parser.add_argument("--max-rb", dest="maxrb", type=int, deault=500, help="The maximum number of reduced bases that will be produced [default: %(default)s]")
parser.add_argument("--num-cores", dest="numcores", type=int, default=1, help="The number of CPUs to request and use [default: %(default)s]")
parser.add_argument("--request-mem", dest="requestmem", type=int, default=1024, help="The amount of RAM (Mb) to request [default: %(default)s]")
parser.add_argument("--max-enrich", dest="maxenrich", type=int, default=1, help="The maximum number of enrichment steps to try [default: %(default)s]")

# parse input options
opts = parser.parse_args()

# make run directory if required
if not os.path.isdir(opts.rundir):
  try:
    os.makedirs(opts.rundir)
  except:
    print("Error... could not make run directory '{}'".format(opts.rundir), file=sys.stderr)
    sys.exit(1)

# check executables
if not os.path.isdir(opts.execdir):
  print("Error... executables directory '{}' does not exist".format(opts.execdir), file=sys.stderr)
  sys.exit(1)
else:
  # set greedycpp executable
  greedyexec = os.path.join(opts.execdir, 'greedyOMP')
  if not os.path.isfile(greedyexec) or not os.access(greedyexec, os.X_OK):
    print("Error... '{}' does not exist or is not executable".format(greedyexec), file=sys.stderr)
    sys.exit(1)

  # set EIM empirical interpolant exectable
  eimexec = os.path.join(opts.execdir, 'eim')
  if not os.path.isfile(eimexec) or not os.access(eimexec, os.X_OK):
    print("Error... '{}' does not exist or is not executable".format(eimexec), file=sys.stderr)
    sys.exit(1)

  # set verfiy 'validation' exectable
  verifyexec = os.path.join(opts.execdir, 'verifyOMP')
  if not os.path.isfile(verifyexec) or not os.access(verifyexec, os.X_OK):
    print("Error... '{}' does not exist or is not executable".format(verifyexec), file=sys.stderr)
    sys.exit(1)


# create initial training set of RA and declination values (uniform across the sky)
tsfile = os.path.join(opts.rundir, "TS_points.txt")
if opts.ntraining < 1:
  print("Error... number of training points must be greater than zero", file=sys.stderr)
  sys.exit(1)

ra = 2.*np.pi*np.random.rand(1,opts.ntraining)
dec = -(np.pi/2.) + np.arccos(2.*np.random.rand(1,opts.ntraining) - 1.)
try:
  np.savetxt(tsfile, np.concatenate((ra, dec)).T)
except:
  print("Error... could not output training points file '{}'".format(tsfile), file=sys.stderr)
  sys.exit(1)

# set configuration file
cfgfile = os.path.join(rundir, "run.cfg")
cfgdic = {}
cfgdic["tsfile"] = tsfile
cfgdic["outdir"] = rundir
if opts.starttime > 0:
  cfgdic["starttime"] = opts.starttime
else:
  print("Error... training set start time '{}' is a negative value".format(opts.starttime), file=sys.stderr)
  sys.exit(1)
if opts.endtime > 0 and opts.endtime > opts.starttime:
  cfgdic["endtime"] = opts.endtime
else:
  print("Error... training set end time '{}' is not valid".format(opts.endtime), file=sys.stderr)
  sys.exit(1)
if opts.npoints > 0:
  cfgdic["npoints"] = opts.npoints
else:
  print("Error... training set length '{}' is not valid".format(opts.npoints), file=sys.stderr)
  sys.exit(1)
if opts.detector is in ['H1', 'L1', 'H2', 'V1', 'G1']:
  cfgdic["detector"] = opts.detector
else:
  print("Error... detector '{}' is not valid".format(opts.detector), file=sys.stderr)
  sys.exit(1)
if opts.ephem is in ['DE405', 'DE200', 'DE414', 'DE421']:
  cfgdic["ephemeris"] = opts.ephem
else:
  print("Error... ephemeris '{}' is not valid".format(opts.ephem), file=sys.stderr)
  sys.exit(1)
if opts.units is in ['TCB', 'TDB']:
  cfgdic["timeunits"] = opts.units
else:
  print("Error... time units '{}' is not valid".format(opts.units), file=sys.stderr)
  sys.exit(1)
if opts.tolerance > 0. and opts.tolerance < 1.:
  cfgdic["tolerance"] = opts.tolerance
else:
  print("Error... tolerance '{}' is not valid".format(opts.tolerance), file=sys.stderr)
  sys.exit(1)
if opts.maxrb > 0:
  cfgdic["maxrb"] = opts.maxrb
else:
  print("Error... maximum number of reduced bases '{}' is not valid".format(opts.maxrb), file=sys.stderr)
  sys.exit(1)

try:
  fp = open(cfgfile, "w")
  fp.write(configtemplate.format(**cfgdic))
  fp.close()
except:
  print("Error... problem writing out configuration file '{}'".format(cfgfile), file=sys.stderr)
  sys.exit(1)


# create submit file for greedycpp
greedysub = os.path.join(rundir, "greedy.sub")
greedytemplate = """universe = vanilla
executable = {greedyexec}
request_memory = {reqmem}
machine_count = {ncpus}
arguments = "{greedyargs}"
accounting_group_user = matthew.pitkin
accounting_group = aluk.dev.o3.cw.targeted.rom
log = {log}
error = {error}
output = {output}
getenv = true
notifications = Never
queue 1
"""
greedydic = {}
greedydic["greedyexec"] = greedyexec
if opts.requestmem > 0:
  greedydic["reqmem"] = opts.requestmem
else:
  print("Error... invalid amount of memory '{}' requested".format(opts.requestmem), file=sys.stderr)
  sys.exit(1)
if opts.numcores > 0:
  greedydic["ncpus"] = opts.numcores
else:
  print("Error... invalid number of CPUs '{}' requested".format(opts.numcores), file=sys.stderr)
  sys.exit(1)
greedydic["greedyargs"] = cfgfile

logdir = os.path.join(rundir, "log")
if not os.path.isdir(logdir):
  try:
    os.makedirs(logdir)
  except:
    print("Error... could not create log directory '{}'".format(logdir), file=sys.stderr)
    sys.exit(1)
greedydic["log"] = os.path.join(logdir, "greedy.log")
greedydic["error"] = os.path.join(logdir, "greedy.err")
greedydic["output"] = os.path.join(logdir, "greedy.out")

try:
  fpgreedy = open(greedysub, "w")
  fpgreedy.write(greedytemplate.format(**greedydic))
  fpgreedy.close()
except:
  print("Error... could not output greedy '.sub' file '{}'".format(greedysub), file=sys.stderr)
  sys.exit(1)

# create DAG file
dagfile = os.path.join(rundir, "run.dag")
dagfp = open(dagfile, "w")

# create submit file for



