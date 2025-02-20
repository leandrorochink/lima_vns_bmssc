//============================================================================
// Author      : Leandro R. Costa, Daniel Aloise, Nenad Mladenovic.
// Description : Implementation of the LIMA-VNS published in the paper "Less is 
//               more: basic variable neighborhood search heuristic for 
//               balanced minimum sum-of-squares clustering". Please, check
//               https://doi.org/10.1016/j.ins.2017.06.019 for theoretical 
//               details. 
//============================================================================

#include "Vns.h"
#include "Solution.h"
#include "LocalSearch.h"
#include <iostream>
#include <vector>
#include <list>
#include <iterator>
#include <ctime>
#include <cstdlib>
#include <random>
#include "tempsC++.h"
#include <stdlib.h>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "Pair.h"
#include <limits>
#include <algorithm>


using namespace std; 

struct c_unique {
  int current;
  c_unique() {current=0;}
  int operator()() {return current++;}
} UniqueNumber;

Vns::Vns(vector<Point>* _dataset, DistanceMatrix* _distances, int _nClusters, Random* _random, vector< vector<Pair> >* _rankedEntities){
	dataset = _dataset;
	distances = _distances;
	nClusters = _nClusters;
	random = _random;
	rankedEntities = _rankedEntities;
	k = 1;
	timer = ChronoCPU();
}

int Vns::execute(Solution& bestSolution, double tempMax, int kMin, int kStep, int KMax, string outputFileName){
	int iteration = 1;
	int nLocalSearchCalls = 0;
	int kMax = KMax;

	LocalSearch localSearch(dataset, random, rankedEntities);

	initialSolution(bestSolution);
	timer.Reset();
	timer.Start();

	bool upgrade;
	k = kMin;

	while(timer.GetTime() <= tempMax){
		Solution partialSolution(bestSolution);

		upgrade = shaking(partialSolution);
		if(upgrade){
			localSearch.execute(partialSolution, &timer, tempMax, iteration);
			nLocalSearchCalls++;
		}

		if((partialSolution.solutionValue + 0.01) < bestSolution.solutionValue){
			bestSolution.copy(partialSolution);
			k = kMin;
		}else{
			k+=kStep;
		}

		if(k>kMax){
			k = kMin;
			iteration++;
		}
	}
	return iteration;
}

void Vns::initialSolution(Solution& initial){
	vector<int> entities(dataset->size());

	initial.time = 0.0;

	generate(entities.begin(), entities.end(), UniqueNumber);
	random->random_shuffle(entities.begin(), entities.end());

	for(int index=0; index<(signed)dataset->size(); index++){
		int cluster = index%nClusters;
		initial.assignment[entities[index]] = cluster;
		initial.clusterSizes[cluster] += 1;
	}

	initial.solutionValue = 0.0;
	for(unsigned int i=0; i<dataset->size()-1; i++){
		for(unsigned int j=i+1; j<dataset->size(); j++){
			if(initial.assignment[i] == initial.assignment[j]){
				initial.solutionValue += distances->getDistance(i,j)/initial.clusterSizes[initial.assignment[i]];
			}
		}
	}

	initial.initializeSc();
}

bool Vns::shaking(Solution& solution){
	int nData = dataset->size();
	vector <int> choosedEntities(nData, 0);
	int pointI;
	int pointJ;

	for(int l=0; l<k; l++){
		pointI = random->get_rand_ij(0,nData-1);
		pointJ = random->get_rand_ij(0,nData-1);
		while(choosedEntities[pointI] == 1 || choosedEntities[pointJ] == 1 || pointI == pointJ
						|| solution.assignment[pointI] == solution.assignment[pointJ]){
			pointI = random->get_rand_ij(0,nData-1);
			pointJ = random->get_rand_ij(0,nData-1);
		}

		int clusterI = solution.assignment[pointI];
		int clusterJ = solution.assignment[pointJ];
		solution.assignment[pointI] = clusterJ;
		solution.assignment[pointJ] = clusterI;

		solution.solutionValue += solution.sc[pointI][clusterJ]/solution.clusterSizes[clusterJ] - solution.sc[pointI][clusterI]/solution.clusterSizes[clusterI]
							   + solution.sc[pointJ][clusterI]/solution.clusterSizes[clusterI] - solution.sc[pointJ][clusterJ]/solution.clusterSizes[clusterJ]
                               - distances->getDistance(pointI, pointJ)/solution.clusterSizes[clusterI] - distances->getDistance(pointI, pointJ)/solution.clusterSizes[clusterJ];
		for(int i=0; i<nData; i++){
			solution.sc[i][clusterI] += - distances->getDistance(i,pointI) + distances->getDistance(i,pointJ);
			solution.sc[i][clusterJ] += - distances->getDistance(i,pointJ) + distances->getDistance(i,pointI);
		}
	}
	solution.time = timer.GetTime();
	return true;
}

bool Vns::checkSolution(Solution* solution){
	double nClusters = solution->nClusters;
	double nDimentions = dataset->at(0).getDimensions();
	int nData = dataset->size();

	vector<Point> centroids;
	double solutionValueCalculed = 0.0;

	vector< vector<double> > centroidsAux;
	for(int i=0; i<nClusters; i++){
		vector<double> aux(nDimentions, 0.0);
		centroidsAux.push_back(aux);
	}

	for(int i=0; i<nData; i++){
		for(int d=0; d<nDimentions; d++){
			centroidsAux[solution->assignment[i]][d] += (dataset->at(i).getCoordinatesAt(d)/solution->clusterSizes[solution->assignment[i]]);
		}
	}
	vector <int> clustersCount(nClusters, 0);
	for(int i=0; i<nClusters; i++){
		centroids.push_back(Point(centroidsAux[i]));
	}

	for(int i=0; i<nData; i++){
		solutionValueCalculed += centroids[solution->assignment[i]].getSquaredDistance(dataset->at(i));
		clustersCount[solution->assignment[i]] = clustersCount[solution->assignment[i]] + 1;
	}

	double dif = fabs(solution->solutionValue - solutionValueCalculed);
	if(dif >= 0.0001){
		cout << "################################################################################" << endl;
		cout << "DIVERGENCIA NO VALOR DA SOLUCAO!"<<endl;
		cout << "Divergencia = " << setprecision(8) << scientific << dif << endl;
		cout << "ValueSol = " << setprecision(8) << scientific << solution->solutionValue << endl;
		cout << "ValueCal = " << setprecision(8) << scientific << solutionValueCalculed << endl;
		cout << "################################################################################" << endl;
		return false;
	}

	for(int i=0; i<nClusters; i++){
		if(solution->clusterSizes[i] != clustersCount[i]){
			cout << "################################################################################" << endl;
			cout << "DIVERGENCIA NOS TAMANHOS DOS CLUSTERS!"<<endl;
			cout << "Cluster = " << setprecision(2) << fixed << i << endl;
			cout << "################################################################################" << endl;
			return false;
		}
	}
	return true;
}

