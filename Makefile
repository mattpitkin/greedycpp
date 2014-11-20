all: greedy_mpi verify greedy

### model specific flags, headers, sources ###
MODELFLAGS=
MODELSOURCES=
MODELHEADERS=models/spa_waveforms.h
### end model specific information ###

COMPILER=g++
COMPILER_MPI=mpicxx

OPENMP=-fopenmp -DUSE_OPENMP

LDLIBS = -lgsl -lgslcblas -lconfig++
#CXXFLAGS=-g -O0 `gsl-config --cflags`
#LDLIBS = `gsl-config --libs` -L/opt/local/lib -lconfig++
CXXFLAGS=


### gcc optimizations for 64-bit OS on x86_64 with corei7 (see README) ###
#OPTFLAGS=-O3 -ffast-math -ftree-vectorizer-verbose=6 -ftree-vectorize -mavx -DOPTIMIZE_AVX
OPTFLAGS=-O3 -mtune=corei7-avx -ffast-math -ftree-vectorize -DOPTIMIZE_AVX

### icc optimizations (see README) ###
#OPTFLAGS=-O3 -xHOST -DOPTIMIZE_AVX


SOURCES = training_set.cpp parameters.cpp
HEADERS = training_set.hpp parameters.hpp gsl_helper_functions.h  my_models.h utils.h


greedy: greedy.cpp $(SOURCES) $(HEADERS) $(MODELSOURCES) $(MODELHEADERS)
	$(COMPILER) $(CXXFLAGS) $(OPTFLAGS) $(MODELFLAGS) $(OPENMP) -DCOMPILE_WITHOUT_MPI -o greedy greedy.cpp $(SOURCES) $(MODELSOURCES) $(LDLIBS)

greedy_mpi: greedy.cpp $(SOURCES) $(HEADERS) $(MODELSOURCES) $(MODELHEADERS)
	$(COMPILER_MPI) $(CXXFLAGS) $(OPTFLAGS) $(MODELFLAGS) $(OPENMP) -o greedympi greedy.cpp $(SOURCES) $(MODELSOURCES) $(LDLIBS)

verify: basis_validation.cpp $(SOURCES) $(HEADERS) $(MODELSOURCES) $(MODELHEADERS)
	$(COMPILER) $(CXXFLAGS) $(MODELFLAGS)  $(OPTFLAGS) $(OPENMP) -o verify basis_validation.cpp $(SOURCES) $(MODELSOURCES) $(LDLIBS)

.PHONY: clean
clean:
	rm -f greedy greedympi verify


