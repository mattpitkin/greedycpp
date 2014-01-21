all: GreedyCpp

GreedyCpp: greedy.cpp
	g++ -o greedy greedy.cpp -lgsl -lgslcblas
