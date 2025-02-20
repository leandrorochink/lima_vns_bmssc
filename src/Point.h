//============================================================================
// Author      : Leandro R. Costa, Daniel Aloise, Nenad Mladenovic.
// Description : Implementation of the LIMA-VNS published in the paper "Less is 
//               more: basic variable neighborhood search heuristic for 
//               balanced minimum sum-of-squares clustering". Please, check
//               https://doi.org/10.1016/j.ins.2017.06.019 for theoretical 
//               details. 
//============================================================================

#ifndef POINT_H_
#define POINT_H_

#include <string>
#include <vector>

using namespace std;

class Point{
    vector<double> coordinates;

public:
	Point(){}
    Point(vector<double> _coordinates);
	double getCoordinatesAt(int index);
	void setCoordinatesAt(int index, double value);
	const vector<double> getCoordinates();
    double getDistance(Point point);
    double getSquaredDistance(Point point);
    double getSquaredDistance(vector<double> coord);
	int getDimensions();

    string toString();
};
#endif
