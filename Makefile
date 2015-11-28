SHELL=/bin/bash

### record git hash ###
GIT_SHA1=$(shell git rev-parse HEAD)
$(info "the git hash is"  ${GIT_SHA1})
GIT_FLAGS=-DGIT_SHA1=\"$(GIT_SHA1)\"

### set your c++ compiler(s) here ###
CXX=g++
CXX_MPI=mpicxx

### source and executable directories shouldn't change ###
SRCDIR=code
BINDIR=bin

### OpenMP flags ###
OPENMP=-fopenmp -DUSE_OPENMP

### Work with numpy data files ###
NUMPY=-DUSE_NUMPY
NUMPYHEADERS=-I/home/balzani57/pool/include
NUMPYLIBS=-L/home/balzani57/pool/lib -lcnpy

### GSL and GSL's verion of blas  ###
#LDLIBS = -lgsl -lgslcblas -lconfig++ -pg
LDLIBS = `gsl-config --libs` -L/opt/local/lib
CXXFLAGS=`gsl-config --cflags`
#CXXFLAGS=-g -O0 `gsl-config --cflags` $(shell pkg-config --cflags gsl && pkg-config --cflags lalsimulation)

### Input file parser libconfig++ ###
LDLIBS+=-lconfig++
CXXFLAGS+=
#CXXFLAGS=-g -O0 `gsl-config --cflags` $(shell pkg-config --cflags gsl && pkg-config --cflags lalsimulation)


### Optimizations ###
# gcc optimizations for 64-bit OS on x86_64 with corei7 (see README) #
#OPTFLAGS=-O3 -ffast-math -ftree-vectorizer-verbose=6 -ftree-vectorize -mavx -DOPTIMIZE_AVX
OPTFLAGS=-O3 -mtune=corei7-avx -ffast-math -ftree-vectorize -DOPTIMIZE_AVX
# icc optimizations (see README) #
#OPTFLAGS=-O3 -xHOST -DOPTIMIZE_AVX


### model specific flags, headers, sources ###
MODELFLAGS=
MODELSOURCES=
MODELHEADERS=models/taylorf2/spa_waveforms.h
MODELLIBS=

## LIGO Analysis Library (LAL) flags for LAL models ##
LAL=  # if you do not have LAL installed
#LAL=-DMODEL_LAL 
#MODELHEADERS+=models/lal/phenomp.h
#MODELFLAGS+=$(shell pkg-config --cflags lalsimulation)
#MODELLIBS+=$(shell pkg-config --libs lalsimulation)




### Source and headers ###
SOURCES = $(SRCDIR)/training_set.cpp \
	$(SRCDIR)/parameters.cpp \
	$(SRCDIR)/gsl_helper_functions.cpp \
	$(SRCDIR)/quadratures.cpp \
	$(SRCDIR)/load_simulation_data.cpp \
	$(SRCDIR)/eim.cpp

HEADERS = $(SRCDIR)/training_set.hpp \
	$(SRCDIR)/parameters.hpp \
	$(SRCDIR)/gsl_helper_functions.hpp \
	$(SRCDIR)/quadratures.hpp \
	$(SRCDIR)/load_simulation_data.hpp \
	$(SRCDIR)/eim.hpp \
	$(SRCDIR)/my_models.h \
	$(SRCDIR)/utils.h

all: $(BINDIR)/greedy_mpi \
	$(BINDIR)/greedyOMP_mpi \
	$(BINDIR)/verifyOMP \
	$(BINDIR)/greedy \
	$(BINDIR)/greedyOMP \
	$(BINDIR)/eim

$(BINDIR)/greedy: $(SRCDIR)/greedy.cpp $(SOURCES) $(HEADERS) $(MODELSOURCES) $(MODELHEADERS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(MODELFLAGS) $(LAL) $(NUMPY) $(GIT_FLAGS) \
        $(NUMPYHEADERS) -DCOMPILE_WITHOUT_MPI -o $(BINDIR)/greedy $(SRCDIR)/greedy.cpp \
        $(SOURCES) $(MODELSOURCES) $(LDLIBS) $(MODELLIBS) $(NUMPYLIBS)

$(BINDIR)/greedyOMP: $(SRCDIR)/greedy.cpp $(SOURCES) $(HEADERS) $(MODELSOURCES) $(MODELHEADERS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(MODELFLAGS) $(OPENMP) $(LAL) $(NUMPY) $(GIT_FLAGS) \
        $(NUMPYHEADERS) -DCOMPILE_WITHOUT_MPI -o $(BINDIR)/greedyOMP $(SRCDIR)/greedy.cpp \
	$(SOURCES) $(MODELSOURCES) $(LDLIBS) $(MODELLIBS) $(NUMPYLIBS)

$(BINDIR)/greedy_mpi: $(SRCDIR)/greedy.cpp $(SOURCES) $(HEADERS) $(MODELSOURCES) $(MODELHEADERS)
	$(CXX_MPI) $(CXXFLAGS) $(OPTFLAGS) $(MODELFLAGS) $(LAL) $(NUMPY) $(GIT_FLAGS) \
        $(NUMPYHEADERS) -o $(BINDIR)/greedy_mpi $(SRCDIR)/greedy.cpp $(SOURCES) \
	$(MODELSOURCES) $(LDLIBS) $(MODELLIBS) $(NUMPYLIBS)

$(BINDIR)/greedyOMP_mpi: $(SRCDIR)/greedy.cpp $(SOURCES) $(HEADERS) $(MODELSOURCES) $(MODELHEADERS)
	$(CXX_MPI) $(CXXFLAGS) $(OPTFLAGS) $(MODELFLAGS) $(OPENMP) $(LAL) $(NUMPY) $(GIT_FLAGS) \
        $(NUMPYHEADERS) -o $(BINDIR)/greedyOMP_mpi $(SRCDIR)/greedy.cpp $(SOURCES) \
	$(MODELSOURCES) $(LDLIBS) $(MODELLIBS) $(NUMPYLIBS)

$(BINDIR)/verifyOMP: $(SRCDIR)/validation.cpp $(SOURCES) $(HEADERS) $(MODELSOURCES) $(MODELHEADERS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(MODELFLAGS) $(OPENMP) $(LAL) $(NUMPY) \
	$(NUMPYHEADERS) -o $(BINDIR)/verifyOMP $(SRCDIR)/validation.cpp $(SOURCES) \
	$(MODELSOURCES) $(LDLIBS) $(MODELLIBS) $(NUMPYLIBS)

$(BINDIR)/eim: $(SRCDIR)/eim_driver.cpp $(SOURCES) $(HEADERS) $(MODELSOURCES) $(MODELHEADERS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(MODELFLAGS) $(OPENMP) $(LAL) $(NUMPY) \
	$(NUMPYHEADERS) -o $(BINDIR)/eim $(SRCDIR)/eim_driver.cpp $(SOURCES) \
	$(MODELSOURCES) $(LDLIBS) $(MODELLIBS) $(NUMPYLIBS)

.PHONY: clean
clean:
	rm -f bin/greedy bin/greedyOMP bin/greedy_mpi bin/greedyOMP_mpi bin/verifyOMP bin/eim

test: test1 test2 test3 test4

test1:
	mkdir -p Temporary
	./bin/greedy examples/example1/test1.cfg > Temporary/Test1.txt
	./scripts/regression_test.py MyOutputTest1/ examples/output_test1/
	rm -rf MyOutputTest1/
test2:
	mkdir -p Temporary
	mpirun -np 3 bin/greedy_mpi examples/example1/test1.cfg > Temporary/Test2.txt
	./scripts/regression_test.py MyOutputTest1/ examples/output_test1/
	rm -rf MyOutputTest1/
test3:
	mkdir -p Temporary
	./bin/greedyOMP examples/example1/test1.cfg > Temporary/Test3.txt
	./scripts/regression_test.py MyOutputTest1/ examples/output_test1/
	rm -rf MyOutputTest1/
test4:
	mkdir -p Temporary
	mpirun -np 1 bin/greedyOMP_mpi examples/example1/test1.cfg > Temporary/Test4.txt
	./scripts/regression_test.py MyOutputTest1/ examples/output_test1/
	rm -rf MyOutputTest1/test: test1 test2 test3 test4

test1:
	mkdir -p Temporary
	./bin/greedy examples/example1/test1.cfg > Temporary/Test1.txt
	./scripts/regression_test.py MyOutputTest1/ examples/output_test1/
	rm -rf MyOutputTest1/
test2:
	mkdir -p Temporary
	mpirun -np 3 bin/greedy_mpi examples/example1/test1.cfg > Temporary/Test2.txt
	./scripts/regression_test.py MyOutputTest1/ examples/output_test1/
	rm -rf MyOutputTest1/
test3:
	mkdir -p Temporary
	./bin/greedyOMP examples/example1/test1.cfg > Temporary/Test3.txt
	./scripts/regression_test.py MyOutputTest1/ examples/output_test1/
	rm -rf MyOutputTest1/
test4:
	mkdir -p Temporary
	mpirun -np 1 bin/greedyOMP_mpi examples/example1/test1.cfg > Temporary/Test4.txt
	./scripts/regression_test.py MyOutputTest1/ examples/output_test1/
	rm -rf MyOutputTest1/test: test1 test2 test3 test4

test1:
	mkdir -p Temporary
	./bin/greedy examples/example1/test1.cfg > Temporary/Test1.txt
	./scripts/regression_test.py MyOutputTest1/ examples/output_test1/
	rm -rf MyOutputTest1/
test2:
	mkdir -p Temporary
	mpirun -np 3 bin/greedy_mpi examples/example1/test1.cfg > Temporary/Test2.txt
	./scripts/regression_test.py MyOutputTest1/ examples/output_test1/
	rm -rf MyOutputTest1/
test3:
	mkdir -p Temporary
	./bin/greedyOMP examples/example1/test1.cfg > Temporary/Test3.txt
	./scripts/regression_test.py MyOutputTest1/ examples/output_test1/
	rm -rf MyOutputTest1/
test4:
	mkdir -p Temporary
	mpirun -np 1 bin/greedyOMP_mpi examples/example1/test1.cfg > Temporary/Test4.txt
	./scripts/regression_test.py MyOutputTest1/ examples/output_test1/
	rm -rf MyOutputTest1/
