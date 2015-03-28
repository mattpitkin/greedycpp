SHELL=/bin/bash

### set your c++ compiler(s) here ###
CXX=g++
CXX_MPI=mpicxx

### source and executable directories shouldn't change ###
SRCDIR=code
BINDIR=bin

### model specific flags, headers, sources ###
MODELFLAGS=$(shell pkg-config --cflags lalsimulation)
MODELSOURCES=
MODELHEADERS=models/taylorf2/spa_waveforms.h models/lal/phenomp.h
MODELLIBS=$(shell pkg-config --libs lalsimulation)

### OpenMP flags ###
OPENMP=-fopenmp -DUSE_OPENMP

### Needed for gsl, gsl's verion of blas and input file parser libconfig++ ###
#LDLIBS = -lgsl -lgslcblas -lconfig++ -pg
LDLIBS = `gsl-config --libs` -L/opt/local/lib -lconfig++
CXXFLAGS=
#CXXFLAGS=-g -O0 `gsl-config --cflags` $(shell pkg-config --cflags gsl && pkg-config --cflags lalsimulation)

### Optimizations ###
# gcc optimizations for 64-bit OS on x86_64 with corei7 (see README) #
#OPTFLAGS=-O3 -ffast-math -ftree-vectorizer-verbose=6 -ftree-vectorize -mavx -DOPTIMIZE_AVX
OPTFLAGS=-O3 -mtune=corei7-avx -ffast-math -ftree-vectorize -DOPTIMIZE_AVX
# icc optimizations (see README) #
#OPTFLAGS=-O3 -xHOST -DOPTIMIZE_AVX

### Source and headers ###
SOURCES = $(SRCDIR)/training_set.cpp $(SRCDIR)/parameters.cpp
HEADERS = $(SRCDIR)/training_set.hpp $(SRCDIR)/parameters.hpp \
          $(SRCDIR)/gsl_helper_functions.h $(SRCDIR)/my_models.h $(SRCDIR)/utils.h


all: $(BINDIR)/greedy_mpi $(BINDIR)/greedyOMP_mpi $(BINDIR)/verifyOMP $(BINDIR)/greedy $(BINDIR)/greedyOMP

$(BINDIR)/greedy: $(SRCDIR)/greedy.cpp $(SOURCES) $(HEADERS) $(MODELSOURCES) $(MODELHEADERS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(MODELFLAGS) \
        -DCOMPILE_WITHOUT_MPI -o $(BINDIR)/greedy $(SRCDIR)/greedy.cpp \
        $(SOURCES) $(MODELSOURCES) $(LDLIBS) $(MODELLIBS)

$(BINDIR)/greedyOMP: $(SRCDIR)/greedy.cpp $(SOURCES) $(HEADERS) $(MODELSOURCES) $(MODELHEADERS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(MODELFLAGS) $(OPENMP) \
        -DCOMPILE_WITHOUT_MPI -o $(BINDIR)/greedyOMP $(SRCDIR)/greedy.cpp \
	$(SOURCES) $(MODELSOURCES) $(LDLIBS) $(MODELLIBS)

$(BINDIR)/greedy_mpi: $(SRCDIR)/greedy.cpp $(SOURCES) $(HEADERS) $(MODELSOURCES) $(MODELHEADERS)
	$(CXX_MPI) $(CXXFLAGS) $(OPTFLAGS) $(MODELFLAGS) \
        -o $(BINDIR)/greedy_mpi $(SRCDIR)/greedy.cpp $(SOURCES) \
	$(MODELSOURCES) $(LDLIBS) $(MODELLIBS)

$(BINDIR)/greedyOMP_mpi: $(SRCDIR)/greedy.cpp $(SOURCES) $(HEADERS) $(MODELSOURCES) $(MODELHEADERS)
	$(CXX_MPI) $(CXXFLAGS) $(OPTFLAGS) $(MODELFLAGS) $(OPENMP) \
        -o $(BINDIR)/greedyOMP_mpi $(SRCDIR)/greedy.cpp $(SOURCES) \
	$(MODELSOURCES) $(LDLIBS) $(MODELLIBS)

$(BINDIR)/verifyOMP: $(SRCDIR)/basis_validation.cpp $(SOURCES) $(HEADERS) $(MODELSOURCES) $(MODELHEADERS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(MODELFLAGS) $(OPENMP) \
        -o $(BINDIR)/verifyOMP $(SRCDIR)/basis_validation.cpp $(SOURCES) \
	$(MODELSOURCES) $(LDLIBS) $(MODELLIBS)

.PHONY: clean
clean:
	rm -f bin/greedy bin/greedyOMP bin/greedy_mpi bin/greedyOMP_mpi bin/verifyOMP


