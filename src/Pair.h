//============================================================================
// Author      : Leandro R. Costa, Daniel Aloise, Nenad Mladenovic.
// Description : Implementation of the LIMA-VNS published in the paper "Less is 
//               more: basic variable neighborhood search heuristic for 
//               balanced minimum sum-of-squares clustering". Please, check
//               https://doi.org/10.1016/j.ins.2017.06.019 for theoretical 
//               details. 
//============================================================================

#ifndef PAIR_H_
#define PAIR_H_

class Pair {
public:
	Pair(int _id, double _value);

	double getValue() const;
	void setValue(double _value);
	int getId() const;
	void setId(int _id);

	bool operator<(const Pair& other) const{
		if(this->value < other.getValue()){
			return true;
		}
		return false;
	};

private:
	int id;
	double value;
};

#endif /* PAIR_H_ */
