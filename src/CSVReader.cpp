//============================================================================
// Author      : Leandro R. Costa, Daniel Aloise, Nenad Mladenovic.
// Description : Implementation of the LIMA-VNS published in the paper "Less is 
//               more: basic variable neighborhood search heuristic for 
//               balanced minimum sum-of-squares clustering". Please, check
//               https://doi.org/10.1016/j.ins.2017.06.019 for theoretical 
//               details. 
//============================================================================
#include "CSVReader.h"
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include "Point.h"
#include <vector>

vector<Point> Reader::readInstance(string pathFile){
	std::ifstream file(pathFile.c_str());

	std::vector<std::vector<double> > ints;

	bool insertNewLine = true;

	double oneInt;

	while (file.good() && file >> oneInt){
		if(file.eof()){
			break;
		}

		if (insertNewLine){
			std::vector<double> vc;
			ints.push_back(vc);

			insertNewLine = false;
		}

		ints.back().push_back(oneInt);

		int ch;

		while ((ch = file.peek()) != std::char_traits<char>::eof()){
			if (ch == ','|| ch == '\t' || ch == '\r' || ch == '\n'){
				char ch2;

				if (!file.read(&ch2, 1)){
					break;
				}

				if (ch == '\n' && !insertNewLine){
					insertNewLine = true;
				}
			}else{
				break;
			}
		}
	}

	vector<Point> dataset;

	for(int i=0; i<(signed)ints.size(); i++){
		vector<double> coordinates;
		for(int j=0; j<(signed)ints[i].size(); j++){
			coordinates.push_back(ints[i][j]);
		}
		dataset.push_back(coordinates);
	}
	return dataset;
}

vector< vector<double> > Reader::readTimesFile(string pathFile){
	vector< vector<double> > times(16);
	vector<double> aux(10, 0.0);

	for(int i=0; i<16; i++){
		times[i] = aux;
	}

	char delimiter = ';';
	ifstream file(pathFile.c_str());

	string value;
	for(int i=0; i<16; i++ ){
		for(int j=0; j<9; j++){
			getline (file, value, delimiter);

			if(file.eof() && !file.good()){
				break;
			}
			times[i][j] = atof(value.c_str());
		}
		getline (file, value, '\n');
		if(file.eof() && !file.good()){
			break;
		}
		times[i][9] = atof(value.c_str());
	}
	return times;
}

const char* Reader::returnPrintable(string value){
	string aux = "";
	for (unsigned int i = 0; i < value.length(); i++) {
		if (isprint(value[i])) {
			aux += value[i];
		}
	}
	return aux.c_str();
}
