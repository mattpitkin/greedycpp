all: greedy_mpi verify

LDLIBS = -lgsl -lgslcblas -lconfig++ -lhdf5

### comment out COMPILE_WITH_MPI defined in greedy.cpp ###
greedy: greedy.cpp
	g++ -o greedy greedy.cpp training_set.cpp parameters.cpp $(LDLIBS)

greedy_mpi: greedy.cpp
	mpicxx -o greedympi greedy.cpp training_set.cpp parameters.cpp $(LDLIBS)

verify: basis_validation.cpp
	g++ -o verify basis_validation.cpp parameters.cpp training_set.cpp -lconfig++ -lgsl -lgslcblas

.PHONY: clean
clean:
	rm -f greedy greedympi verify


