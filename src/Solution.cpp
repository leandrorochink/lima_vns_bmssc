//============================================================================
// Author      : Leandro R. Costa, Daniel Aloise, Nenad Mladenovic.
// Description : Implementation of the LIMA-VNS published in the paper "Less is 
//               more: basic variable neighborhood search heuristic for 
//               balanced minimum sum-of-squares clustering". Please, check
//               https://doi.org/10.1016/j.ins.2017.06.019 for theoretical 
//               details. 
//============================================================================

#include "Solution.h"
#include "DistanceMatrix.h"
#include <vector>
#include <iostream>

using namespace std;

Solution::Solution(){
	 solutionValue = 0;
 }

Solution::Solution(int _nClusters, int _nDataPoints, DistanceMatrix* _distances){
	nClusters = _nClusters;
	nDataPoints = _nDataPoints;
	distances =  _distances;
	solutionValue = 0;
	time = 0.0;

	assignment = new int[nDataPoints];
	clusterSizes = new double[nClusters];

	sc = new double*[nDataPoints];
	for(int i = 0; i<nDataPoints; i++){
		sc[i] = new double[nClusters];
	}

	for(int i=0; i<nClusters; i++){
		clusterSizes[i] = 0;
	}

}

Solution::Solution(const Solution& copy){
	distances = copy.distances;
	nClusters = copy.nClusters;
	nDataPoints = copy.nDataPoints;
	time = copy.time;
	solutionValue = copy.solutionValue;

	assignment = new int[nDataPoints];
	clusterSizes = new double[nClusters];

	sc = new double*[nDataPoints];
	for(int i=0; i<nDataPoints; i++){
		assignment[i] = copy.assignment[i];
		sc[i] = new double[nClusters];
		for(int j=0; j<nClusters; j++){
			sc[i][j] = copy.sc[i][j];
		}
	}

	for(int i=0; i<nClusters; i++){
		clusterSizes[i]=copy.clusterSizes[i];
	}
}

Solution::~Solution(){
	for(int i=0; i<nDataPoints; i++){
		delete [] sc[i];
	}
	delete [] sc;
	delete [] assignment;
	delete [] clusterSizes;
}

void Solution::copy(const Solution& copy){
	distances = copy.distances;
	nClusters = copy.nClusters;
	nDataPoints = copy.nDataPoints;
	time = copy.time;
	solutionValue = copy.solutionValue;

	for(int i=0; i<nDataPoints; i++){
		assignment[i] = copy.assignment[i];
		for(int j=0; j<nClusters; j++){
			sc[i][j] = copy.sc[i][j];
		}
	}

	for(int i=0; i<nClusters; i++){
		clusterSizes[i]=copy.clusterSizes[i];
	}
}

void Solution::initializeSc(){
	for(int i=0; i<nDataPoints; i++){
		for(int j=0; j<nClusters; j++){
			sc[i][j] = 0.0;
		}
	}

	for(int i=0; i<nDataPoints; i++){
		for(int j=0; j<nDataPoints; j++){
			sc[i][assignment[j]] += distances->getDistance(i,j);
		}
	}
}














