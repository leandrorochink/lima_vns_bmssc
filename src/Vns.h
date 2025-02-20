//============================================================================
// Author      : Leandro R. Costa, Daniel Aloise, Nenad Mladenovic.
// Description : Implementation of the LIMA-VNS published in the paper "Less is 
//               more: basic variable neighborhood search heuristic for 
//               balanced minimum sum-of-squares clustering". Please, check
//               https://doi.org/10.1016/j.ins.2017.06.019 for theoretical 
//               details. 
//============================================================================

#ifndef VNS_H_
#define VNS_H_

#include "Solution.h"
#include "Point.h"
#include "DistanceMatrix.h"
#include "LocalSearch.h"
#include <vector>
#include <random>
#include "Random.h"
#include <iomanip>
#include "Pair.h"

using namespace std;

class Vns {
public:
	Vns(vector<Point>* _dataset, DistanceMatrix* _distances, int _nClusters, Random* _random, vector< vector<Pair> >* _rankedEntities);
	int execute(Solution& bestSolution, double tempMax, int kMin, int kStep, int kMax, string outputFileName);

private:
	int nClusters;
	int k;

	Random* random;
	vector<Point>* dataset;
	DistanceMatrix* distances;
	vector< vector<Pair> >* rankedEntities;

	ChronoCPU timer;

	bool shaking(Solution& solution);
	void initialSolution(Solution& initial);
	bool checkSolution(Solution* solution);
};
#endif /* VNS_H_ */
