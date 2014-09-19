all: greedy_mpi

LDLIBS = -lgsl -lgslcblas -lconfig++ -lhdf5

### comment out COMPILE_WITH_MPI defined in greedy.cpp ###
greedy: greedy.cpp
	g++ -o greedy greedy.cpp training_set.cpp parameters.cpp training_space.cpp $(LDLIBS)

greedy_mpi: greedy.cpp
	mpicxx -o greedympi greedy.cpp training_set.cpp parameters.cpp training_space.cpp $(LDLIBS)

.PHONY: clean
clean:
	rm -f greedy greedympi


