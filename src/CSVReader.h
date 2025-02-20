//============================================================================
// Author      : Leandro R. Costa, Daniel Aloise, Nenad Mladenovic.
// Description : Implementation of the LIMA-VNS published in the paper "Less is 
//               more: basic variable neighborhood search heuristic for 
//               balanced minimum sum-of-squares clustering". Please, check
//               https://doi.org/10.1016/j.ins.2017.06.019 for theoretical 
//               details. 
//=============================================================================

#ifndef CSVREADER_H
#define	CSVREADER_H

#include <string>
#include "Point.h"

class Reader {
public:
	vector<Point> readInstance(std::string file);
	vector< vector<double> > readTimesFile(string pathFile);
private:
    const char* returnPrintable(string value);
};
#endif	/* CSVREADER_H */
