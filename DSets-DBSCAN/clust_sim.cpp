//
//  clust_sim.cpp
//  DSets-DBSCAN
//
//  Created by Rui Ma on 29.06.18.
//  Copyright © 2018 Rui Ma. All rights reserved.
//
//this function is for generate simulated 2D cluster data
//the simulated cluster is a spherical cluster with a constant noise
//inputs:
//1) int n_clust: number of clusters
//2) int n_point: number of points in each cluster
//3) int n_noise: number of noise
//4) double r: radius of cluster
//4) double height: number of rows of the image matrix
//5) double width: number of column of the image matrix
//output:
//1) vector<double> clust_data:
//                  col1: x pos
//                  col2: y pos
//                  col3: cluster ID, 0 indicates noise
///////////////////////////////////////////////////////
#define _USE_MATH_DEFINES

#include "clust_sim.hpp"
#include <iostream>
#include <math.h>
#include <armadillo>

using namespace std;
using namespace arma;

mat clust_sim(int n_clust, int n_point, int n_noise, double r, double height, double width){
    //declare output
    mat clust_data(n_point * n_clust + n_noise, 2);
    //generate clusters
    for (int i = 0; i < n_clust; i++){
		clust_data(span(i * n_point, (i + 1) * n_point), span(0, 1)) 
			= sim_spherical_clust(n_point, r);
		clust_data(span(i * n_point, (i + 1) * n_point), span(2, 2))
			= i;
    }
	for (int i = n_point * n_clust; i < n_point * n_clust + n_noise; i++) {
		clust_data(i, 0) = height * randu();
		clust_data(i, 1) = width * randu();
		clust_data(i, 1)
			= 0;
	}
    return clust_data;
}

mat sim_spherical_clust(int n_point, double r){
    //col: x y
    mat clust_data(n_point, 2);
    mat theta = (2 * M_PI) * randu(n_point, 1);
    mat phi = acos(1 - 2 * randu(n_point, 1));
    clust_data.col(0) = sin(phi) % cos(theta);
    clust_data.col(1) = sin(phi) % sin(theta);
    for (int i = 0; i < n_point; i++){
        clust_data(i, 0) = r * clust_data(i, 0) + 0.1 * r * (randu() - 0.5);
        clust_data(i, 1) = r * clust_data(i, 1) + 0.1 * r * (randu() - 0.5);
    }
    return clust_data;
}
