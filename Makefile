all: GreedyMPI
##Greedy GreedyMPI
#
Greedy: greedy.cpp
	g++ -o greedy greedy.cpp TrainingSet.cpp -lgsl -lgslcblas

GreedyMPI: greedy.cpp
	mpicxx -o greedympi greedy.cpp TrainingSet.cpp -lgsl -lgslcblas -lconfig++
