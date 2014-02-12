all: GreedyMPI
##Greedy GreedyMPI
#
Greedy: greedy.cpp
	g++ -o greedy greedy.cpp -lgsl -lgslcblas

GreedyMPI: greedy.cpp
	mpicxx -o greedympi greedy.cpp -lgsl -lgslcblas
