//============================================================================
// Author      : Leandro R. Costa, Daniel Aloise, Nenad Mladenovic.
// Description : Implementation of the LIMA-VNS published in the paper "Less is 
//               more: basic variable neighborhood search heuristic for 
//               balanced minimum sum-of-squares clustering". Please, check
//               https://doi.org/10.1016/j.ins.2017.06.019 for theoretical 
//               details. 
//============================================================================

#include "LocalSearch.h"
#include "Solution.h"
#include <vector>
#include "tempsC++.h"
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <stdlib.h>
#include "Pair.h"
#include <map>
#include <random>

using namespace std;

LocalSearch::LocalSearch(vector<Point>* _dataset, Random* _random, vector< vector<Pair> >* _rankedEntities){
	dataset = _dataset;
	random = _random;
	rankedEntities = _rankedEntities;
}


void LocalSearch::execute(Solution& bestLocalSolution, ChronoCPU* timer, double maxTime, int nIteration){
	bool upgrade = true;
	while(upgrade){
		upgrade = false;
		upgrade = swapLocalSearchFirstRand(bestLocalSolution, timer, maxTime);
	}
}

bool LocalSearch::swapLocalSearchBest(Solution& solution, ChronoCPU* timer, double maxTime){
	bool upgrade = false;
	int bestI = 0;
	int bestJ = 0;
	int bestCI = 0;
	int bestCJ = 0;
	double bestDf = 0.0;
	double bestTime = 0.0;

	for(unsigned int i=0; i<dataset->size()-1; i++){
		for(unsigned int j=i+1; j<dataset->size(); j++){
			if(solution.assignment[i] != solution.assignment[j]){

				int clusterI = solution.assignment[i];
				int clusterJ = solution.assignment[j];

				double df = (solution.sc[i][clusterJ] - solution.sc[j][clusterJ] - solution.distances->getDistance(i,j))/solution.clusterSizes[clusterJ]
						   + (solution.sc[j][clusterI] - solution.sc[i][clusterI] - solution.distances->getDistance(i,j))/solution.clusterSizes[clusterI];

				double time = timer->GetTime();
				if(time > maxTime){
					i = dataset->size()-1;
					break;
				}

				if(df < bestDf){
					upgrade = true;
					bestI = i;
					bestJ = j;
					bestCI = clusterI;
					bestCJ = clusterJ;
					bestDf = df;
					bestTime = time;



				}
			}
		}
	}
	if(upgrade){
		swap(solution, bestCI, bestI, bestCJ, bestJ, bestDf);
		solution.time = bestTime;
		return true;
	}else{
		return false;
	}
}

bool LocalSearch::swapLocalSearchFirstRand(Solution& solution, ChronoCPU* timer, double maxTime){
	int start = random->get_rand_ij(0, dataset->size()-1);
	for(unsigned int i = start; i < dataset->size()+start-1; i++){
		for(unsigned int j=i+1; j < dataset->size()+start; j++){

			if(solution.assignment[i%dataset->size()] != solution.assignment[j%dataset->size()]){
				int clusterI = solution.assignment[i%dataset->size()];
				int clusterJ = solution.assignment[j%dataset->size()];

				double df = (solution.sc[i%dataset->size()][clusterJ] - solution.sc[j%dataset->size()][clusterJ] - solution.distances->getDistance(i%dataset->size(),j%dataset->size()))/solution.clusterSizes[clusterJ]
                          + (solution.sc[j%dataset->size()][clusterI] - solution.sc[i%dataset->size()][clusterI] - solution.distances->getDistance(i%dataset->size(),j%dataset->size()))/solution.clusterSizes[clusterI];


				double time = timer->GetTime();
				if(time > maxTime){
					return false;
				}

				if(df < 0.0){
					swap(solution, clusterI, i%dataset->size(), clusterJ, j%dataset->size(), df);
					solution.time = time;
					return true;
				}
			}
		}
	}
	return false;
}

void LocalSearch::swap(Solution& solution, int clusterI, int i, int clusterJ, int j, double df){

	solution.assignment[i] = clusterJ;
	solution.assignment[j] = clusterI;

	solution.solutionValue += df;

	for(int k=0; k<solution.nDataPoints; k++){
		solution.sc[k][clusterI] += - solution.distances->getDistance(k,i) + solution.distances->getDistance(k,j);
		solution.sc[k][clusterJ] += - solution.distances->getDistance(k,j) + solution.distances->getDistance(k,i);
	}
}


bool LocalSearch::checkSolution(Solution* solutionBefore, Solution* solutionAfter, double deltaSolutionValue){
	double nClusters = solutionBefore->nClusters;
	double nDimentions = dataset->at(0).getDimensions();
	int nData = dataset->size();

	double solutionValueCalculedBefore = 0.0;
	double solutionValueCalculedAfter = 0.0;

	vector<int> clustersCountB(nClusters, 0);
	vector<int> clustersCountA(nClusters, 0);

	vector< vector<double> > centroidsAuxB;
	vector< vector<double> > centroidsAuxA;
	for(int i=0; i<nClusters; i++){
		vector<double> aux(nDimentions, 0.0);
		centroidsAuxB.push_back(aux);
		centroidsAuxA.push_back(aux);
	}

	for(int i=0; i<nData; i++){
		for(int d=0; d<nDimentions; d++){
			centroidsAuxB[solutionBefore->assignment[i]][d] += (dataset->at(i).getCoordinatesAt(d)/solutionBefore->clusterSizes[solutionBefore->assignment[i]]);
			centroidsAuxA[solutionAfter->assignment[i]][d] += (dataset->at(i).getCoordinatesAt(d)/solutionAfter->clusterSizes[solutionAfter->assignment[i]]);
		}
	}

	for(int i=0; i<nData; i++){
		solutionValueCalculedBefore += dataset->at(i).getSquaredDistance(centroidsAuxB[solutionBefore->assignment[i]]);
		solutionValueCalculedAfter += dataset->at(i).getSquaredDistance(centroidsAuxA[solutionAfter->assignment[i]]);
		clustersCountB[solutionBefore->assignment[i]] = clustersCountB[solutionBefore->assignment[i]] + 1;
		clustersCountA[solutionAfter->assignment[i]] = clustersCountA[solutionAfter->assignment[i]] + 1;
	}

	double difB = fabs(solutionBefore->solutionValue - solutionValueCalculedBefore);
	if(difB >= 0.0001){
		cout << "################################################################################" << endl;
		cout << "DIVERGENCIA NO VALOR DA SOLUCAO ANTIGA!"<<endl;
		cout << "Divergencia = " << setprecision(8) << scientific << difB << endl;
		cout << "ValueSol = " << setprecision(8) << scientific << solutionBefore->solutionValue << endl;
		cout << "ValueCal = " << setprecision(8) << scientific << solutionValueCalculedBefore << endl;
		cout << "################################################################################" << endl;
		return false;
	}

	double difA = fabs(solutionAfter->solutionValue - solutionValueCalculedAfter);
	if(difA >= 0.0001){
		cout << "################################################################################" << endl;
		cout << "DIVERGENCIA NO VALOR DA SOLUCAO NOVA!"<<endl;
		cout << "Divergencia = " << setprecision(8) << scientific << difA << endl;
		cout << "ValueSol = " << setprecision(8) << scientific << solutionAfter->solutionValue << endl;
		cout << "ValueCal = " << setprecision(8) << scientific << solutionValueCalculedAfter << endl;
		cout << "################################################################################" << endl;
		return false;
	}

	double dfCalculated = fabs(solutionValueCalculedBefore - solutionValueCalculedAfter);
	double dif = fabs(dfCalculated - fabs(deltaSolutionValue));
	if(dif >= 0.0001){
		cout << "################################################################################" << endl;
		cout << "DIVERGENCIA NO VALOR DO DELTA!"<<endl;
		cout << "Divergencia = " << setprecision(8) << scientific << dif << endl;
		cout << "DeltaLS = " << setprecision(8) << scientific << deltaSolutionValue << endl;
		cout << "DeltaCalc = " << setprecision(8) << scientific << dfCalculated << endl;
		cout << "################################################################################" << endl;
		return false;
	}

	for(int i=0; i<nClusters; i++){
		if(solutionBefore->clusterSizes[i] != clustersCountB[i]){
			cout << "################################################################################" << endl;
			cout << "DIVERGENCIA NOS TAMANHOS DOS CLUSTERS DA SOLUCAO ANTIGA!"<<endl;
			cout << "Cluster = " << setprecision(2) << fixed << i << endl;
			cout << "################################################################################" << endl;
			return false;
		}
	}

	for(int i=0; i<nClusters; i++){
		if(solutionAfter->clusterSizes[i] != clustersCountA[i]){
			cout << "################################################################################" << endl;
			cout << "DIVERGENCIA NOS TAMANHOS DOS CLUSTERS DA SOLUCAO NOVA!"<<endl;
			cout << "Cluster = " << setprecision(2) << fixed << i << endl;
			cout << "################################################################################" << endl;
			return false;
		}
	}
	return true;
}
