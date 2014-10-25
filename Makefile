all: greedy_mpi verify

LDLIBS = -lgsl -lgslcblas -lconfig++ -lhdf5
#CXXFLAGS=-g -O0 `gsl-config --cflags`
CXXFLAGS=
#LDLIBS = -lgsl -lgslcblas -lconfig++ -lhdf5
#LDLIBS = `gsl-config --libs` -L/opt/local/lib -lconfig++ -lhdf5

SOURCES = training_set.cpp parameters.cpp
HEADERS = training_set.hpp parameters.hpp gauss_wgts.h spa_waveforms.h gsl_helper_functions.h  my_models.h utils.h

### comment out COMPILE_WITH_MPI defined in greedy.cpp ###
greedy: $(SOURCES) $(HEADERS)
	g++ $(CXXFLAGS) -o greedy greedy.cpp $(SOURCES) $(LDLIBS)

greedy_mpi: greedy.cpp $(SOURCES) $(HEADERS)
	mpicxx $(CXXFLAGS) -o greedympi greedy.cpp $(SOURCES) $(LDLIBS)

verify: basis_validation.cpp $(SOURCES) $(HEADERS)
	g++ $(CXXFLAGS) -o verify basis_validation.cpp $(SOURCES) -lconfig++ -lgsl -lgslcblas

.PHONY: clean
clean:
	rm -f greedy greedympi verify


