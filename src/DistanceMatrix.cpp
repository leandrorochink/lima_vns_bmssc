//============================================================================
// Author      : Leandro R. Costa, Daniel Aloise, Nenad Mladenovic.
// Description : Implementation of the LIMA-VNS published in the paper "Less is 
//               more: basic variable neighborhood search heuristic for 
//               balanced minimum sum-of-squares clustering". Please, check
//               https://doi.org/10.1016/j.ins.2017.06.019 for theoretical 
//               details. 
//============================================================================

#include "DistanceMatrix.h"
#include <iostream>
#include <vector>
#include "Point.h"

DistanceMatrix::DistanceMatrix(vector<Point>* dataset){
	nV = dataset->size();

	adj = new double*[nV];
	for(int i=0; i<nV; i++){
		adj[i] = new double[nV-i];
	}

	for(int i=0; i<nV; i++){
		setDistance(i, i, 0.0);
		for(int j=i+1; j<nV; j++){
			setDistance(i, j, (*dataset)[i].getSquaredDistance((*dataset)[j]));
		}
	}
}

DistanceMatrix::~DistanceMatrix(){
	for(int i=0; i<nV; i++){
		delete [] adj[i];
	}
	delete [] adj;
}

double DistanceMatrix::getDistance(int i, int j){
	if(i<j){
		return adj[i][j-i];
	}else{
		return adj[j][i-j];
	}
}

void DistanceMatrix::setDistance(int i, int j, double d){
	if(i<j){
		adj[i][j-i] = d;
	}else{
		adj[j][i-j] = d;
	}
}
