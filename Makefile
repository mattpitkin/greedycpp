all: GreedyMPI

LDLIBS = -lgsl -lgslcblas -lconfig++ -lhdf5

### this will not work unless you edit the code -- see README ###
Greedy: greedy.cpp
	g++ -o greedy greedy.cpp training_set.cpp $(LDLIBS)

GreedyMPI: greedy.cpp
	mpicxx -o greedympi greedy.cpp training_set.cpp training_space.cpp $(LDLIBS)

.PHONY: clean
clean:
	rm -f greedy greedympi


