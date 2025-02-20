//============================================================================
// Author      : Leandro R. Costa, Daniel Aloise, Nenad Mladenovic.
// Description : Implementation of the LIMA-VNS published in the paper "Less is 
//               more: basic variable neighborhood search heuristic for 
//               balanced minimum sum-of-squares clustering". Please, check
//               https://doi.org/10.1016/j.ins.2017.06.019 for theoretical 
//               details. 
//============================================================================

#ifndef DISTANCEMATRIX_H
#define DISTANCEMATRIX_H

#include <vector>
#include "Point.h"

using namespace std;

class DistanceMatrix{
    int nV;
    double **adj;

public:
    DistanceMatrix(vector<Point>* dataset);
	~DistanceMatrix();
    double getDistance(int i, int j);
    void setDistance(int i, int j, double d);
};
#endif

