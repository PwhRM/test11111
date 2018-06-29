//
//  DSets_DBSCAN.cpp
//  DSets-DBSCAN
//
//  Created by Rui Ma on 29.06.18.
//  Copyright © 2018 Rui Ma. All rights reserved.
//
//
//DSets-DBSCAN method for clustering
//original papers:
//[1] DSets-DBSCAN: A Parameter-Free Clustering Algorithm
//[2] Infection and immunization:  a new class of evolutionarygame dynamics.
//
//inputs:
//1) mat dist_mat: Euclidean distance map, arma mat matrix
//2) mat x0: initial population state. If omottted or empty then x0 will be located in the baricenter of the n-dimensional symplex
//3) double supportThreshold: the threshold used to extract the support from the population vector x (default 1e-4)
//4) double precision: the maximum population distance between two successive steps to consider the dynamics in equilibrium (default 1e-8)
//5) maxIters: maximum number of iterations of the dynamic systems (default 1000)
//output:
//1) mat clust: arma mat matrix, clust id
///////////////////////////////////////////////////////
#include "DSets_DBSCAN.hpp"
#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

mat DSets_DBSCAN(){
    mat clust;
    
    return clust;
}
