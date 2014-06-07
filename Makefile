all: GreedyMPI

LDLIBS = -lgsl -lgslcblas -lconfig++ -lhdf5

Greedy: greedy.cpp
	g++ -o greedy greedy.cpp TrainingSet.cpp $(LDLIBS)

GreedyMPI: greedy.cpp
	mpicxx -o greedympi greedy.cpp TrainingSet.cpp $(LDLIBS)
