all: GreedyMPI

LDLIBS = -lgsl -lgslcblas -lconfig++ -lhdf5

### this will not work unless you edit the code -- see README ###
Greedy: greedy.cpp
	g++ -o greedy greedy.cpp TrainingSet.cpp $(LDLIBS)

GreedyMPI: greedy.cpp
	mpicxx -o greedympi greedy.cpp TrainingSet.cpp $(LDLIBS)
