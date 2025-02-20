//============================================================================
// Author      : Leandro R. Costa, Daniel Aloise, Nenad Mladenovic.
// Description : Implementation of the LIMA-VNS published in the paper "Less is 
//               more: basic variable neighborhood search heuristic for 
//               balanced minimum sum-of-squares clustering". Please, check
//               https://doi.org/10.1016/j.ins.2017.06.019 for theoretical 
//               details. 
//============================================================================

#include "Pair.h"

using namespace std;

Pair::Pair(int _id, double _value){
	id = _id ;
	value = _value;
}

double Pair::getValue() const {
	return value;
}

void Pair::setValue(double _value) {
	this->value = _value;
}

int Pair::getId() const {
	return id;
}

void Pair::setId(int _id) {
	this->id = _id;
}


