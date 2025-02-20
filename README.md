# LIMA-VNS for Balanced Minimum Sum-of-Squares Clustering
Datasets and C++ source code of the LIMA-VNS algorithm published in the paper "Less is more: basic variable neighborhood search heuristic for balanced minimum sum-of-squares clustering" by Leandro R. Costa, Daniel Aloise, Nenad Mladenovic. 

Please, check https://doi.org/10.1016/j.ins.2017.06.019 for theoretical details. 

# Compiling

You only have to use the makfile indire the *"src"* folder:

*make*

To clean the compiled files use:

*make clean*

# Executing

32-bits version:

*./lima_vns_32 <path/instance.csv> <k=number of clusters> <cpu time limit> <number of runs> <seed> <path/output file> <path/cluster assignment file>*

64-bits version:

*./lima_vns_64 <path/instance.csv> <k=number of clusters> <cpu time limit> <number of runs> <seed> <path/output file> <path/cluster assignment file>*

The cluster assignments file has the instance used as the first column. The remaining columns are the cluster assignments, for example, "iris.csv, 0, 0, 1, ...", which means: 
 - point 1 is assigned to cluster 0; 
 - point 2 is assigned to cluster 0; 
 - point 3 is assigned to cluster 1;

