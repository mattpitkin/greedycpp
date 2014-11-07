all: greedy_mpi verify

LDLIBS = -lgsl -lgslcblas -lconfig++
#CXXFLAGS=-g -O0 `gsl-config --cflags`
CXXFLAGS=
#LDLIBS = -lgsl -lgslcblas -lconfig++
#LDLIBS = `gsl-config --libs` -L/opt/local/lib -lconfig++

SOURCES = training_set.cpp parameters.cpp
HEADERS = training_set.hpp parameters.hpp gauss_wgts.h gsl_helper_functions.h  my_models.h utils.h

### model specific flags, headers, sources ###
MODELFLAGS=
MODELSOURCES=
MODELHEADERS=models/spa_waveforms.h

### comment out COMPILE_WITH_MPI defined in greedy.cpp ###
greedy: $(SOURCES) $(HEADERS) $(MODELSOURCES) $(MODELHEADERS)
	g++ $(CXXFLAGS) $(MODELFLAGS) -o greedy greedy.cpp $(SOURCES) $(MODELSOURCES) $(LDLIBS)

greedy_mpi: greedy.cpp $(SOURCES) $(HEADERS) $(MODELSOURCES) $(MODELHEADERS)
	mpicxx $(CXXFLAGS) $(MODELFLAGS) -o greedympi greedy.cpp $(SOURCES) $(MODELSOURCES) $(LDLIBS)

verify: basis_validation.cpp $(SOURCES) $(HEADERS) $(MODELSOURCES) $(MODELHEADERS)
	g++ $(CXXFLAGS) $(MODELFLAGS) -o verify basis_validation.cpp $(SOURCES) $(MODELSOURCES) -lconfig++ -lgsl -lgslcblas

.PHONY: clean
clean:
	rm -f greedy greedympi verify


