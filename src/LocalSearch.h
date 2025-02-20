//============================================================================
// Author      : Leandro R. Costa, Daniel Aloise, Nenad Mladenovic.
// Description : Implementation of the LIMA-VNS published in the paper "Less is 
//               more: basic variable neighborhood search heuristic for 
//               balanced minimum sum-of-squares clustering". Please, check
//               https://doi.org/10.1016/j.ins.2017.06.019 for theoretical 
//               details. 
//============================================================================

#ifndef LOCALSEARCH_H_
#define LOCALSEARCH_H_
#include "Solution.h"
#include "tempsC++.h"
#include "Pair.h"
#include <random>
#include "Random.h"
#include "Pair.h"

using namespace std;

class LocalSearch {
private:
	vector<Point>* dataset;
	Random* random;
	vector< vector<Pair> >* rankedEntities;

public:

	LocalSearch(vector<Point>* _dataset, Random* _random, vector< vector<Pair> >* _rankedEntities);
	void execute(Solution& bestLocalSolution, ChronoCPU* timer, double maxTime, int nIteration);
	bool swapLocalSearchBest(Solution& solution, ChronoCPU* timer, double maxTime);
	bool swapLocalSearchFirstRand(Solution& solution, ChronoCPU* timer, double maxTime);
	void swap(Solution& solution, int i, int pointV, int j, int pointU, double df);

	bool checkSolution(Solution* solutionBefore, Solution* solutionAfter, double deltaSolutionValue);
};
#endif /* LOCALSEARCH_H_ */
