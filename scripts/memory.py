"""
Module to estimate memory usage 
""" 
    
import os, re, sys
import numpy as np


def estimate_memory(ts_size,param_dim,max_rb,num_workers,quad_size):

    sizeof_double = 8.0e-9 # double in GB

    memory_worker = sizeof_double*( ts_size*param_dim + num_workers + 2*max_rb*ts_size/num_workers + 2*quad_size*ts_size/num_workers )

    memory_master = 2*sizeof_double*( num_workers + quad_size*max_rb + max_rb*max_rb) + sizeof_double*ts_size*param_dim

    #print 'Each worker processor requires at least %1.2f GB' % memory_worker
    print 'All worker processor require at least %1.2f GB (in total)' % (memory_worker*num_workers)
    print 'Master (orthogonalization) processor requires at least %1.2f GB' % memory_master
    print 'worker memory estimate very accurate. Master about factor of 5 too small. Why?'


if __name__=="__main__":

    if ( len(sys.argv) != 6):
        raise Exception("Must specify ts_size, param_dim, max_RB, numprocs_worker, xQuad_size")

    ts_size     = float(sys.argv[1])
    param_dim   = float(sys.argv[2])
    max_rb      = float(sys.argv[3])
    num_workers = float(sys.argv[4])
    quad_size   = float(sys.argv[5])

    estimate_memory(ts_size,param_dim,max_rb,num_workers,quad_size)

