//============================================================================
// Author      : Leandro R. Costa, Daniel Aloise, Nenad Mladenovic.
// Description : Implementation of the LIMA-VNS published in the paper "Less is 
//               more: basic variable neighborhood search heuristic for 
//               balanced minimum sum-of-squares clustering". Please, check
//               https://doi.org/10.1016/j.ins.2017.06.019 for theoretical 
//               details. 
//============================================================================

#include <iostream>
#include "Vns.h"
#include "DistanceMatrix.h"
#include <string>
#include "CSVReader.h"
#include <float.h>
#include <iomanip>
#include <fstream>
#include "Random.h"
#include <sstream>
#include "Pair.h"
#include <algorithm>

using namespace std;

int main(int argc, char** argv) {
	Reader reader;
	vector<Point> dataset;
	double bestSolutionValue = DBL_MAX;
	double bestTime;

	////////////// Parameters //////////////

	string path_instance;
	int n_clusters;
	double max_time;
	int n_runs;
	int seed;
	string path_output;
	string path_output_assignment;

	///////////////////////////////////////

	if(argc < 8){
		cout << "ARGUMENT(S) MISSING!!" << endl << "Usage: " << argv[0];
		cout << " <path/instance.csv> <k=number of clusters> <cpu time limit>";
		cout << " <number of runs> <seed> <path/output file> <path/assignment file>" << endl;
		return EXIT_FAILURE;
	}else{
		 path_instance = argv[1];
		 n_clusters = atoi(argv[2]);
		 max_time = atof(argv[3]);
		 n_runs = atoi(argv[4]);
		 seed = atoi(argv[5]);
		 path_output = argv[6];
		 path_output_assignment = argv[7];
	}



	try{
		ifstream instance_file(path_instance.c_str(), std::ifstream::in);
		if(!instance_file.good()){
			cout << "PROBLEM IN THE PATH OF THE INSTANCE FILE" << endl;
			return EXIT_FAILURE;
		}
	}catch(...){
			cerr << "PROBLEM IN THE PATH OF THE INSTANCE FILE" << endl;
			return EXIT_FAILURE;
	}


	stringstream str_stats;
	str_stats << path_output << ".csv";
	try{
		ofstream results_stats_file(str_stats.str().c_str(), ofstream::app);
		if(!results_stats_file.good()){
			cout << "PROBLEM IN THE PATH OF THE OUTPUT FILE" << endl;
			return EXIT_FAILURE;
		}
	}catch(...){
		cerr << "PROBLEM IN THE PATH OF THE OUTPUT FILE" << endl;
		return EXIT_FAILURE;
	}
	ofstream results_stats_file;
	results_stats_file.open(str_stats.str().c_str(), ofstream::app);



	stringstream str_assignment;
	str_assignment << path_output_assignment << ".csv";

	try{
		ofstream results_assignment_file(str_assignment.str().c_str(), ofstream::app);

		if(!results_assignment_file.good()){
			cout << "PROBLEM IN THE PATH OF THE ASSIGNMENT FILE" << endl;
			return EXIT_FAILURE;
		}
	}catch(...){
		cerr << "PROBLEM IN THE PATH OF THE ASSIGNMENT FILE" << endl;
		return EXIT_FAILURE;
	}
	ofstream results_assignment_file;
	results_assignment_file.open(str_assignment.str().c_str(), ofstream::app);


    int kMin = 2;
	int averageVnsIteration = 0;
	dataset = reader.readInstance(path_instance);
	DistanceMatrix distances(&dataset);
	Solution bestSolution(n_clusters, dataset.size(), &distances);

	vector< vector<Pair> > rankedEntities(dataset.size());
	for(unsigned int o=0; o<dataset.size(); o++){
		for(unsigned int m=0; m<dataset.size(); m++){
			if(o!=m){
				Pair pair(m, distances.getDistance(o,m));
				rankedEntities[o].push_back(pair);
			}
		}
		sort(rankedEntities[o].begin(), rankedEntities[o].end());
	}

	int kMax = dataset.size()/2;
	int kStep = (int)kMax/20;

	cout << "============================================================================================================" << endl;
	cout << "Instance: " << path_instance << endl;
	cout << "Clusters: " << n_clusters << endl;
	cout << "Kmax: " << kMax << endl;
	cout << "KStep: " << kStep<< endl;

	double mean = 0.0;
	double timeMean = 0.0;
	bestSolutionValue = DBL_MAX;
	for(int j=0; j<n_runs; j++){
		Random random(seed);
		Vns vns(&dataset, &distances, n_clusters, &random, &rankedEntities);

		cout << "------------------------------------- Execution " <<  j+1 << " -----------------------------------------" << endl;
		cout << "Seed = " << seed << endl;
		cout << "maxTime = " << setprecision(4) << fixed << max_time << endl;
		Solution solution(n_clusters, dataset.size(), &distances);
		stringstream ss;

		int nIteration = vns.execute(solution, max_time, kMin, kStep, kMax, ss.str());

		if(solution.solutionValue < bestSolutionValue){
			bestSolution.copy(solution);
			bestSolutionValue = solution.solutionValue;
			bestTime = solution.time;
		}
		averageVnsIteration += nIteration;
		timeMean += solution.time;
		mean += solution.solutionValue;

		cout << endl << setprecision(8)<< scientific << "Objective Function value: ";
		cout << solution.solutionValue << " in " << setprecision(4) << fixed << solution.time;
		cout << " seconds"<< endl;

		seed += 1;
	}
	cout <<endl<<"**************************************************************************************"<<endl<<endl;
	cout << "Best Objective Function value found: " << setprecision(8) << scientific << bestSolutionValue << " in " << setprecision(4) << fixed<< bestTime << " seconds"<< endl;
	cout << "Average Objective Function value: " << setprecision(8)<<scientific<< mean/n_runs << endl;
	cout << "Average Time value: "<< setprecision(4) <<fixed<< timeMean/n_runs <<"s"<< endl << endl;

	results_stats_file << path_instance;
	results_stats_file << "," << setprecision(8) << scientific << bestSolutionValue;
	results_stats_file << "," << setprecision(8) << scientific << mean/n_runs;
	results_stats_file << "," << setprecision(4) << fixed << bestTime;
	results_stats_file << "," << setprecision(4) << fixed << timeMean/n_runs << endl;

	results_assignment_file << path_instance;
	for(int i=0; i < (signed)dataset.size(); i++){
		results_assignment_file << ','  << bestSolution.assignment[i];
	}
	results_assignment_file << endl;

	cout <<"**************************************************************************************"<<endl;

	results_stats_file.close();
	results_assignment_file.close();
	return 0;
}





























