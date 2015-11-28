""" Module to estimate memory usage.

The following heap memory allocations comprise the majority of memory usage.
Listed are those that scale like number of processors, max_rb, and ts_size*param_dim.

Format... (setting) file: allocation

Allocated on each worker proc
-----------------------------

GreedyCpp/training_set.cpp:   params_    = new double*[ts_size_];
GreedyCpp/training_set.cpp:   params_[j] = new double[param_dim_];
size(params_) = sizeof(double)*ts_size*param_dim

GreedyCpp/training_set.cpp:    mystart_         = new int[numprocs_worker];
GreedyCpp/training_set.cpp:    myend_           = new int[numprocs_worker];
GreedyCpp/training_set.cpp:    matrix_sub_size_ = new int[numprocs_worker];
size(above 3) = 3*sizeof(int)*numprocs_worker ~ sizeof(double)*numprocs_worker

GreedyCpp/greedy.cpp:    project_coeff = gsl_matrix_complex_alloc(max_RB,ts->matrix_sub_size()[rank]);
size(projec_coeff) = 2*sizeof(double)*max_RB*matrix_sub_size()[rank-1] 

GreedyCpp/greedy.cpp:    TS_gsl = gsl_matrix_complex_alloc(ptspace_class->ts_size(),xQuad->size);
GreedyCpp/greedy.cpp:    TS_gsl = gsl_matrix_complex_alloc(ptspace_class->matrix_sub_size()[rank-1],xQuad->size);
Size(TS_gsl) = 2*sizeof(double)*xQuad_size*matrix_sub_size()[rank-1] 

Total on worker proc = sizeof(double)*[ts_size*param_dim + numprocs_worker + 2*max_RB*matrix_sub_size()[rank-1] + 2*xQuad_size*matrix_sub_size()[rank-1] ]

Since matrix_sub_size()[rank-1] is about ts_size/numprocs_worker

Total on worker proc = sizeof(double)*[ts_size*param_dim + numprocs_worker + 2*max_RB*ts_size/numprocs_worker + 2*xQuad_size*ts_size/numprocs_worker ]

and whenever ts_size >> numprocs_worker we have 

Total on worker proc = sizeof(double)*ts_size[param_dim + 2*max_RB/numprocs_worker + 2*xQuad_size/numprocs_worker ]


Allocated on master proc
------------------------

GreedyCpp/greedy.cpp:    worst_workers_mpi = new int[size];
GreedyCpp/greedy.cpp:    worst_errs_mpi    = new double[size];
size(above 2) ~ 2*sizeof(double)*(numprocs_worker)

GreedyCpp/greedy.cpp:    RB_space          = gsl_matrix_complex_alloc(max_RB,cols);
Size(RB_space) = 2*sizeof(double)*xQuad_size*max_RB 

GreedyCpp/greedy.cpp:    R_matrix          = gsl_matrix_complex_alloc(max_RB,max_RB);
Size(R_matrix) = 2*sizeof(double)*max_RB^2

AND ts class from worker proc section is also built on head node

Total on master proc = 2*sizeof(double)*[numprocs_worker + xQuad_size*max_RB + max_RB^2] + sizeof(double)*ts_size*param_dim



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

