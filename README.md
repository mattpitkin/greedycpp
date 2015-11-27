Overview
--------

A fast, scalable and easy-to-use parallelization of the greedy algorithm
for building an application-specific reduced basis. 

In many cases, for example when implemented on a computer and using
the Euclidean inner product, the greedy algorithm is QR with column pivoting.
And so this code may also be used to compute QR decompositions of a large
matrix.

For details on the code's scaling and QR-based model reduction:

[1] Harbir Antil, Dangxing Chen, and Scott Field. "QR-based model 
    reduction for large problems". In preparation. 


Code's coordinates: https://sfield83@bitbucket.org/sfield83/greedycpp.git


Dependencies
------------

Requirements: libconfig, gsl and MPI (see below to build with mpi disabled).
Optional:     cnpy for reading and writing to Python's numpy binary format.

The get-libconfig and get-cnpy scripts (found in the scripts folder) 
will download and install these external packages for you.

cnpy coordinates: https://github.com/rogersce/cnpy.git


Quick start
-----------

         -- serial mode --
>>> make greedy
>>> ./greedy examples/test1.cfg


         -- mpi --
>>> run 'make'
>>> mpirun -np NUM_PROCS ./greedympi examples/test1.cfg


Adding a new model
------------------

Simply make the following changes to my_models.h

1) Add your model-specific header file
2) Invent a descriptive model-specific tag and add
   your model's evaluation (interface routine and tag) to EvaluateModel

Additionally, you might add precompiler directives before and
after your model-specific code so that the code is only
compiled if -DMODEL_NAME is passed to the compiler.


Running example 1
-----------------

Lets run an example problem using the configuration file provided in example1/

>>> mpirun -np 4 ./ greedympi examples/test1.cfg

A directory called 'my_output_test1' contains the output. Check that 
the output is whats expected by doing

>>> python regression_test.py examples/output_test1 my_output_test1/

This should report on the differences as well as generating a summary 
html file. Differences, if they exist, should be machine precision. 

Next, lets modify test1.cfg. The relevant variables should be changed to 

load_from_file = false;
params_num  = [30,30];
params_low  = [1.0,1.0];
params_high = [3.0,3.0];

These settings will, to the first 10-ish digits, reproduce the same training 
set samplings as that give from file "examples/test1_TS_Points.txt". Next, do 

>>> mpirun -np 4 ./ greedympi example1/test1.cfg

You should find minor differences similar to:

>>> python regression_test.py example1/ example1_output/


Mismatches in the following files: 
==>      ApproxErrors.txt : Maximum difference =   9.00391e-14 on line 261
==>      GreedyPoints.txt : Maximum difference =   1.02141e-14 on line 264
==>        Basis_real.txt : Maximum difference =   1.53470e-07 on line 270
==>        Basis_imag.txt : Maximum difference =   1.27563e-07 on line 270
==>            R_real.txt : Maximum difference =   3.39999e-13 on line 270
==>            R_imag.txt : Maximum difference =   4.36700e-13 on line 270


While differences of 1e-07 might seem large, notice that the 269th 
basis is the last one, and its contribution to the approximation 
is on the order of tol = 1e-12


Problem motivation and typical application
------------------------------------------

Given a parameterized model h(x;p), whose physical variable is x and
parameterization p. 

From a collection of model evaluations at known at a discrete training set
of parameter p_i and physical x_j points, we have a training space

   { h(x_j;p_i) }

This training space fills up a matrix A such that the i-th row is the model
evaluation h(\vec{x};p_i). 

The routine greedy.cpp finds a low-rank approximation of A by a 
pivoted QR (aka greedy) decomposition. The output are a collection 
of p_i (greedy parameter values), the orthonormal basis, and an 
R-matrix such that 

   A = R * Q 

where the rows of Q are precisely the basis.


Running on SuperMike2
---------------------

add the following to ~/.soft:
@default
+Intel-13.0.0
@default
+gsl-1.15-Intel-13.0.0
@default
+openmpi-1.6.3-Intel-13.0.0
@default
+fftw-3.3.3-Intel-13.0.

and run 'resoft'. You might need to download and install libconfig++

Using the lsu makefile
1) rename makefiles/Makefile-LSU-Mike2 to "Makefile". run make
2) open scrips/Qsub-lsu and set number of nodes, submit the job


Troubleshooting
---------------

1) If gsl not found run

gsl-config --cflags --libs-without-cblas

and the output says what must be added to g++ to compile and link.

2) On Macs: Install 'libconfig-hr' from macports instead of the package 
   'libconfig' which is just the c-version

3) If ld cannot link to libconfig++ (problem on Supermike2), 
   build locally after getting libconfig

>>> wget http://www.hyperrealm.com/libconfig/libconfig-1.4.9.tar.gz

Then update LD_LIBRARY_PATH to reflect the location of the newly created
libconfig++ shared library.


Optimizations
-------------

The routine WeightedInner comprises the bulk of the computations. We wish to 
optimize this routine with AVX vector instructions. In turn, underlying gsl
functions are called (and gsl_vector_complex_mul is not a BLAS routine), 
and so vectorization would amount to compiling a vectorized version of gsl.

Instead, the following set of tricks (uncovered by Peter Diener) will allow
for vectorization. First the routine WeightedInner is rewritten to directly
work on the gsl data structures. Next, a preprocessor flag OPTIMIZE_AVX 
can be used to decide which version of WeightedInner to use. With OPTIMIZE_AVX
defined, OPTFLAGS (see the makefile) will vectorize this routine.

One should of course check the optimization works for them. Speedup factors
between 1.5 and 2 have been observed in practice. 

Multithreading
--------------

The code can be compiled with OpenMP. Two versions exist: single node/machine
and MPI+openMP. 

ADD MORE HERE!

Large memory jobs
-----------------

For some applications the size of the basis will be large. The basis is
managed by the rank 0 (master) processor which (by default) shares its memory
with worker processors. For large problems it is better to allow the 
master processor to have an entire node. This can be achieved by using a 
hostfile. An example can be found in the scripts folder.

Iterative enrichment
--------------------

Often the training set one wishes to use is too large to process. Or
it is not known ahead of time if the training set constitutes a faithful
sampling of the continuum. Iteratively enriching the basis identifies
pockets of the parameter space where the basis' accuracy is low, and 
(i) combines these "bad" points with previous greedy points then 
(ii) reruns the greedy on the combined set.

The process is partly automated, and can be carried out as 
follows:

1) create a run directory, run_dir
2) create, say, 3 subdirectories called run0, run1 and run2
and another 2 called randset1 and randset2
3) In run0, place all the necessary parameter and data files needed 
for a typical greedy run.
4) Suppose you have access to 10 nodes. After the greedy algorithm
completes you anticipate running 10 separate basis validation studies,
one on each node. Therefore, you should place 10 separate files
into randset1 and 10 separate files into randset2. Each
file should contain a list of verification points (parameter values
at which you will check the basis accuracy).

You are now ready to iteratively refine the basis.

NOTE: make sure the paths set in the configuration
file are absolute! The output_dir location set by the 
configuration file 

1) Submit a greedy job using data in the run0 directory.
2) Verify the basis accuracy using the results written to 
the output directory (whatever you specified in the
output_dir configuration file). This is done using the veryifyOMP 
program -- you would launch 10 separate executables,
each one taking a different random sampling input file
located in randset1
3) A new folder called validations will be created in output_dir,
and there will be 10 separate files containing a list of 
parameter values for which the basis performed poorly.
4) Using the unix cat program, combine these 10 files together
with the GreedyPoints.txt file found in run0. Call this file
TrainingSetIteration1.txt
5) Move TrainingSetIteration1.txt into run1 directory along
with any other files (e.g. the cfg  file) needed to successfully
run the greedy algorithm.
6) Edit the cfg file to use TrainingSetIteration1.txt as the new
training set.

Repeat the above process until your basis has converged (i.e.
it passes the verification step with no bad points).


Scripts
-------

The script generate_trainingpoints.py will generate a trainingpoints file that can be directly used in the greedycpp config file. It requires an input textfile where each column represents a parameter and each row sets a property:
 p1  p2  p3  p4
 min  min  min  min
 max  max  max  max
 dim  dim  dim  dim
 type type type type
An example can be run using the following command from the base directory:
 python scripts/generate_trainingpoints.py examples/params_info_example.txt --out examples/trainingpoints.txt
This will create a file trainingpoints.txt that can be passed on to the config file.
Use --help for more info and options.


Acknowledgements
----------------
Priscilla Canizares, Collin Capano, Michael Pürrer, and Rory Smith for 
careful error reporting, code improvements and testing of early versions of 
the code. Peter Diener for numerous performance optimizations.