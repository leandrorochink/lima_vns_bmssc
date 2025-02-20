//============================================================================
// Author      : Leandro R. Costa, Daniel Aloise, Nenad Mladenovic.
// Description : Implementation of the LIMA-VNS published in the paper "Less is 
//               more: basic variable neighborhood search heuristic for 
//               balanced minimum sum-of-squares clustering". Please, check
//               https://doi.org/10.1016/j.ins.2017.06.019 for theoretical 
//               details. 
//============================================================================

#ifndef SOLUTION_H_
#define SOLUTION_H_

#include <vector>
#include "DistanceMatrix.h"
#include "Point.h"

using namespace std;

class Solution {
public:

	int nClusters;
	int nDataPoints;
	double solutionValue;
	double time;

	DistanceMatrix* distances;

	double** sc;

	int* assignment;
	double* clusterSizes;

	Solution();
	Solution(const Solution& copy);
	Solution(int _nClusters, int _nDataPoints, DistanceMatrix* _distances);
	~Solution();
	void copy(const Solution& copy);
	void initializeSc();
};
#endif /* SOLUTION_H_ */
