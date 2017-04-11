#!/usr/bin/env python
""" From set of greedy points defining an accurate +,x (linear) mode basis,
create a training set of points following the procedure described in Sec 3.E of 

https://arxiv.org/pdf/1604.08253.pdf

This script converts a set of points with a 0,1 polarization tag
(see https://bitbucket.org/sfield83/greedycpp/wiki/LALSimulation%20models.md)
and creates a list of points with a 5 ("hpPLUShcSquared") polarization tag.

Assumes AddReplaceModeDimension.py is located in the same directory as this script.
"""

import subprocess, argparse, shutil, os

def parse_cmd_line():
  "parse command-line arguments"

  parser = argparse.ArgumentParser(description='Setup a quadratic greedy basis directory',\
           formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('--p', type=str,required=True,\
                      help='Absolute path to validated (and accurate) greedy rundir.')

  args=parser.parse_args()
  return args


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if __name__=="__main__":

  args = parse_cmd_line()
  path = args.p

  rundir_quad = path+"/ts_quad"

  if rundir_quad.find('//') !=-1:
    raise ValueError('rundir path cannot contain //')

  if not os.path.exists(rundir_quad):
    os.makedirs(rundir_quad)
  
  subprocess.call("./AddReplaceModeDimension.py --f "+path+"/GreedyPoints.txt -d", shell=True)
  shutil.move(path+"/GreedyPoints_WithoutModes.txt","GreedyPoints_WithoutModes.txt")
  subprocess.call("./AddReplaceModeDimension.py --f GreedyPoints_WithoutModes.txt --m h2", shell=True)
  subprocess.call("mv GreedyPoints_WithoutModes.txt TS_Quad.txt", shell=True)

  print("***Created TS_Quad.txt and run_settings.cfg. Lets run!***")


  new_ts_size = subprocess.check_output(['wc', '-l', 'TS_Quad.txt'])
  new_ts_size = int( new_ts_size.split(' ')[0] )

  old_cfg = path+'/run_settings.cfg'
  fp     = open(old_cfg,'r')
  lines  = fp.readlines()
  fp.close()

  new_cfg = ''

  new_cfg += '// Derived from '+old_cfg+'\n\n\n'


  # Gaurd against multiple occurances of same string (typical when they've been commented out)
  wrote_max_rb  = False
  wrote_output  = False
  wrote_ts_file = False 

  for line in lines:

    if line.find('max_RB') != -1:
      if not wrote_max_rb:
        new_cfg += 'max_RB = '+str(new_ts_size)+';\n'
        wrote_max_rb = True
    elif line.find('output_dir') != -1:
      if not wrote_output:
        new_cfg += 'output_dir = \"'+rundir_quad+'\";\n'
        wrote_output = True
    elif line.find('ts_file') != -1:
      if not wrote_ts_file:
        new_cfg += 'ts_file = \"'+rundir_quad+'/TS_Quad.txt\";\n'
        wrote_ts_file = True
    else:
      new_cfg += line

  fp = open('run_settings.cfg','a')
  fp.write(new_cfg)
  fp.close()

  shutil.move("TS_Quad.txt",rundir_quad)
  shutil.move("run_settings.cfg",rundir_quad)

