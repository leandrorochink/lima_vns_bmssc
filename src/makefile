#============================================================================
# Author      : Leandro R. Costa, Daniel Aloise, Nenad Mladenovic.
# Description : Implementation of the LIMA-VNS published in the paper "Less is 
#               more: basic variable neighborhood search heuristic for 
#               balanced minimum sum-of-squares clustering". Please, check
#               https://doi.org/10.1016/j.ins.2017.06.019 for theoretical 
#               details. 
#============================================================================

CC = g++

TAGS = -Wall -m64 -O3 -std=c++11

OBJS = Point.o Random.o Pair.o CSVReader.o DistanceMatrix.o Solution.o LocalSearch.o Vns.o LIMA_VNS.o

TARGET = lima_vns_64

%.o: %.cpp
	$(CC) $(TAGS) -c -o $@ $< 

all: $(OBJS) 
	$(CC) $(TAGS) -o $(TARGET) $(OBJS) 

clean: 
	rm $(OBJS) $(TARGET) 
