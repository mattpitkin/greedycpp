all: greedy_mpi verify greedy

LDLIBS = -lgsl -lgslcblas -lconfig++
#CXXFLAGS=-g -O0 `gsl-config --cflags`
CXXFLAGS=
#LDLIBS = `gsl-config --libs` -L/opt/local/lib -lconfig++

#OPTFLAGS=-O3 -ffast-math -ftree-vectorizer-verbose=6 -ftree-vectorize -mavx
OPTFLAGS=-O3 -mtune=corei7-avx -ffast-math -ftree-vectorize

AVX_OPT=-DOPTIMIZE_AVX

SOURCES = training_set.cpp parameters.cpp
HEADERS = training_set.hpp parameters.hpp gsl_helper_functions.h  my_models.h utils.h

### model specific flags, headers, sources ###
MODELFLAGS=
MODELSOURCES=
MODELHEADERS=models/spa_waveforms.h

greedy: greedy.cpp $(SOURCES) $(HEADERS) $(MODELSOURCES) $(MODELHEADERS)
	g++ $(CXXFLAGS) $(OPTFLAGS) $(AVX_OPT) $(MODELFLAGS) -DCOMPILE_WITHOUT_MPI -o greedy greedy.cpp $(SOURCES) $(MODELSOURCES) $(LDLIBS)

greedy_mpi: greedy.cpp $(SOURCES) $(HEADERS) $(MODELSOURCES) $(MODELHEADERS)
	mpicxx $(CXXFLAGS) $(OPTFLAGS) $(AVX_OPT) $(MODELFLAGS) -o greedympi greedy.cpp $(SOURCES) $(MODELSOURCES) $(LDLIBS)

verify: basis_validation.cpp $(SOURCES) $(HEADERS) $(MODELSOURCES) $(MODELHEADERS)
	g++ $(CXXFLAGS) $(MODELFLAGS) $(AVX_OPT) -o verify basis_validation.cpp $(SOURCES) $(MODELSOURCES) $(LDLIBS) -fopenmp

.PHONY: clean
clean:
	rm -f greedy greedympi verify


