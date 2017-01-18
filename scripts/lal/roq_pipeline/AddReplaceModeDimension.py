#!/usr/bin/env python

""" For models which support and ``all_parts" tag, this script can
be used to manipulate a text file of parameter values to add,
remove or change the last column entry which describes the
waveform's polarization. 

Tf the model has, say, 6 parameters, the additional 6+1=7th parameter
is for the polarization.

See https://bitbucket.org/sfield83/greedycpp/wiki/LALSimulation%20models.md

This script takes an existing file of points and appends the "mode dimension",
an integer variable taking values of 0,1,2,3,4,5

See the model's header file for the mode <-> ID correspondence. It 
should be: {"+":0,"x":1,"++":2,"xx":3,"+x":4,"h2":5}.

Ex: my_points.txt contains 100 parameter values...  if the desired
modes are "+" and "x", the output file will contain 200 parameter
values -- the first 100 with an extra "0" appended and the last too
with an extra "1" appended."""

import numpy as np
import argparse

def parse_cmd_line():
  "parse command-line arguments"

  parser = argparse.ArgumentParser(description='Append or remove mode(s) dimension to existing training points',\
           formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('--f', type=str,required=True,\
                      help='Path to text file containing traning points.')
  parser.add_argument('--m', type=str,nargs='+',\
                      help='Supply if ADDING model tag type (e.g. + x ++ xx +x h2 for all of them)')
  parser.add_argument('-r', action='store_true',\
                      help='Switch between mode index conventions (see code)')
  parser.add_argument('-d', action='store_true',\
                      help='Delete the mode parameter')

  args=parser.parse_args()
  return args


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if __name__=="__main__":

  np_fmt = '%.14f'
  #np_fmt = '%.6f' # Possibly reason first quad runs failed to produce accurate basis! Retain all digits

  args    = parse_cmd_line()
  infile  = args.f
  Modes   = args.m
  replace = args.r
  delete  = args.d

  if replace:
    assert(Modes is None)
  if Modes is not None:
    assert(not replace)

  parameters_to_modify = np.loadtxt(infile)

  # define the map modes <-> ID
  ModeToID = {"+":0,"x":1,"++":2,"xx":3,"+x":4,"h2":5}

  # define conventions rory <-> scott (add 2)
  RoryToScott = 2

  NewParameterSets = []

  if Modes is not None:
    for mode in Modes:
      print "Adding mode %s with ID %i"%(mode,ModeToID[mode])
      ModeNums = ModeToID[mode] * np.ones((parameters_to_modify.shape[0],1))
      NewParameterSets.append( np.hstack((parameters_to_modify,ModeNums)) )
    ParametersForFile = np.vstack( NewParameterSets )
    #outfile = infile.replace('.txt','_WithModes.txt')
    outfile = infile
  if replace:
    mode_indx     = parameters_to_modify.shape[1]
    without_modes = parameters_to_modify[:,0:(mode_indx-1)]
    modes         = parameters_to_modify[:,-1]
    assert( np.max( np.abs(modes - np.array(modes,dtype=np.int)) ) == 0 ) # mode parameter should be an integer
    for i, mode in enumerate(modes):
      NewParameterSets.append( np.hstack((without_modes[i,:],mode+RoryToScott)) )
    ParametersForFile = np.vstack( NewParameterSets )
    outfile = infile.replace('.txt','_WithoutModes.txt')
    #outfile = infile
  if delete:
    print("removing modes")
    mode_indx     = parameters_to_modify.shape[1]
    without_modes = parameters_to_modify[:,0:(mode_indx-1)]
    ParametersForFile = without_modes
    outfile = infile.replace('.txt','_WithoutModes.txt')

  np.savetxt(outfile,ParametersForFile,fmt=np_fmt)
  print "Created %s with entries..."%(outfile)
  print "shape of new file",ParametersForFile.shape



