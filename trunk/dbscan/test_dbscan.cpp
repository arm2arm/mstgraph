#include "dbscan.h"
#include "distance.h"
#include <iostream>

int main(){

	using namespace Metrics;

	Clustering::Points ps;
	unsigned int ndim;

	// random init points dataset (dims, points)
	//Clustering::randomInit(ps, 2, 100);   
	Clustering::readFromFile(ps,ndim);
	// init: sim threshold, minPts
	Clustering::DBSCAN clusters(ps,0.9, 1); 

	// uniform distribution dataset
	//	clusters.uniformPartition();          

	
	//Distance<Euclidean<Clustering::Point> > de;
	//clusters.computeDistance(de);     
	
	// build similarity  matrix
	Distance<Cosine<Clustering::Point> > d;
	clusters.computeSimilarity(d);     

	// run clustering
	clusters.run_cluster();

	std::cout << clusters;

	int test;
	std::cin >> test;
	return 0;
}