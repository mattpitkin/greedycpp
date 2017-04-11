#!/usr/bin/env python

import subprocess, os, argparse

def parse_cmd_line():
  """parse command-line arguments for uniform grid resample"""

  parser = argparse.ArgumentParser(description=\
            'Setup uniform grid ROQ build from validated, adaptive grid ROQ data',\
           formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('--r', type=str,required=True,\
                      help='output from another greedy run, MUST END WITH / (should contain run_settings.cf and GreedyPoints.txt)')
  parser.add_argument('--min', type=str,required=True,\
                      help='first quadrature node xmin. need to be floats or libcfg throws an exception (ex: 20.0)')
  parser.add_argument('--max', type=str,required=True,\
                      help='last quadrature node xmax. need to be floats or libcfg throws an exception (ex: 4096.0)')
  parser.add_argument('--num', type=int,required=True,\
                      help='number of quadrature points (ex 3 will give [20,2058,4096]')

  args=parser.parse_args()
  return args

if __name__=="__main__":

  args = parse_cmd_line()

  xmin = args.min
  xmax = args.max
  quad_points = args.num

  #rundir = '/panfs/ds08/sxs/sfield/greedy_seob/production_runs/roq_0/'
  rundir = args.r

  assert(len(xmin.split('.'))==2)
  assert(len(xmax.split('.'))==2)

  df = (float(xmax) - float(xmin)) / ( float(quad_points) - 1.0)
  print("df = %f"%df)
  print("1/df = %f"%(1.0/df))

  new_ts = rundir+'GreedyPoints.txt'

  new_ts_size = subprocess.check_output(['wc', '-l', new_ts])
  new_ts_size = int( new_ts_size.split(' ')[0] )

  fp     = open(rundir+'run_settings.cfg','r')
  lines  = fp.readlines()

  fp.close()

  new_cfg = ''

  for line in lines:

    if line.find('quad_type') != -1:
      new_cfg += '//'+line+'\n'
      new_cfg += '// quad settings written by script using 1/df = %f\n'%(1./df)
      new_cfg += 'quad_type = 1;\n'
      new_cfg += 'x_min = '+str(xmin)+';\n'
      new_cfg += 'x_max = '+str(xmax)+';\n'
      new_cfg += 'quad_points = '+str(quad_points)+';\n'
    elif line.find('quad_nodes_file') != -1 or line.find('num_weight_file') != -1:
      new_cfg += '//'+line+'\n'
    elif line.find('max_RB') != -1:
      new_cfg += '//'+line+'\n'
      new_cfg += 'max_RB = '+str(new_ts_size)+';\n'
    elif line.find('ts_file') != -1:
      new_cfg += '//'+line+'\n'
      new_cfg += 'ts_file = \"'+new_ts+'\";\n'
    elif line.find('output_dir') != -1:
      new_cfg += '//'+line+'\n'
      new_cfg += 'output_dir = \"'+os.getcwd()+'\";\n'
    else:
      new_cfg += line

  fp = open('run_settings.cfg','a')
  fp.write(new_cfg)
  fp.close()
