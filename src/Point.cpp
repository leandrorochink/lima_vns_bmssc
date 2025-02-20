//============================================================================
// Author      : Leandro R. Costa, Daniel Aloise, Nenad Mladenovic.
// Description : Implementation of the LIMA-VNS published in the paper "Less is 
//               more: basic variable neighborhood search heuristic for 
//               balanced minimum sum-of-squares clustering". Please, check
//               https://doi.org/10.1016/j.ins.2017.06.019 for theoretical 
//               details. 
//============================================================================

#include "Point.h"
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include <iomanip>

using namespace std;

Point::Point(vector<double> _coordinates){
	coordinates = _coordinates;
}

double Point::getCoordinatesAt(int index){
	return coordinates[index];
}

void Point::setCoordinatesAt(int index, double value){
	coordinates[index] = value;
}

const vector<double> Point::getCoordinates(){
	return coordinates;
}

int Point::getDimensions(){
	return coordinates.size();
}

double Point::getSquaredDistance(Point point){
	double sum = 0.0;
	
	for(unsigned int i=0; i<coordinates.size(); i++){
		//sum += pow((coordinates[i] - point.getCoordinatesAt(i)),2);
		sum += ((coordinates[i] - point.getCoordinatesAt(i))*(coordinates[i] - point.getCoordinatesAt(i)));
	}
	return sum;
}

double Point::getSquaredDistance(vector<double> coord){
	double sum = 0.0;

	for(unsigned int i=0; i<coordinates.size(); i++){
		sum += ((coordinates[i] - coord[i])*(coordinates[i] - coord[i]));
	}
	return sum;
}

double Point::getDistance(Point point){
	double sum = 0.0;
	
	for(unsigned int i=0; i<coordinates.size(); i++){
		//sum += pow((coordinates[i] - point.getCoordinatesAt(i)),2);
		sum += ((coordinates[i] - point.getCoordinatesAt(i))*(coordinates[i] - point.getCoordinatesAt(i)));
	}
	return sqrt(sum);
}

string Point::toString(){
	stringstream str;

	str << "";
	for(unsigned int i=0; i<coordinates.size(); i++){
		str << setprecision(6)<<fixed<< coordinates[i] << " ";
	}

	return str.str();
}
