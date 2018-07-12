//
//  ExpandCluster.cpp
//  DSets-DBSCAN
//
//  Created by Rui Ma on 05.07.18.
//  Copyright © 2018 Rui Ma. All rights reserved.
//
//expand cluster core from DSets method with DBSCAN method
//inputs:
//1) mat dits_mat: arma matrix, Euclidean distance matrix, not the original one but the one changes size in the loop
//2) uword pt: the index of data points in vec d to be visited. This d is a sub-vector of the orginal and would change its size in the while loop. d(i) - 1 is actually the index in the orignal d
//3) int cid: cluster ID for the current cluster
//4) double Eps: Eps used for local dominant
//5) int minPts
//6) vec clust: cluster ID for all the data points, the index of pt is actual d(pt) - 1 in this clust vector
//7) vec d: arma vec, current d in the loop
//output:
//1) mat x: n * 1 arma matrix, the population vector at convergence
///////////////////////////////////////////////////////

#include <iostream>
#include <armadillo>
#include "ExpandCluster.hpp"

using namespace std;
using namespace arma;

vec ExpandCluster(mat dist_mat, uword pt, int cid, double Eps, int minPts, vec clust, vec d){
    //check number of points in the Eps radius but without the corepoint
    //seeds in the row index in dist_mat with col = pt;
    uvec seeds = find(dist_mat.col(pt) <= Eps && dist_mat.col(pt) > 0);
    
    //one exception is the number of points is already smaller than minPts
    if (d.n_elem < minPts - 1) {
        clust(d(pt) - 1) = cid;
        return clust;
    }
    
    //check if noise
    if (seeds.n_elem < minPts - 1) {
        clust(d(pt) - 1) = 0;
        return clust;
    } else {
        //first assign cluster ID
        clust(d(pt) - 1) = cid;
        //core point is already removed from seeds
        //put points in seeds into loop
        while (!seeds.is_empty()) {
            //alway test the first element in the column vector
            //and at the end of the while loop delete the first element
            uword currentP = seeds(0);
            uvec result = find(dist_mat.col(currentP) <= Eps);
            if (result.n_elem >= minPts) {
                for (int i = 0; i < int(result.n_elem); i++) {
                    uword resultP = result(i);
                    if (clust(d(resultP) - 1) == -1 || clust(d(resultP) - 1) == 0) {
                        if (clust(d(resultP) - 1) == -1) {
                            seeds.insert_rows(seeds.n_elem, 1);
                            seeds(seeds.n_elem - 1) = resultP;
                        }
                        //assign cluster id to resultP
                        clust(d(resultP) - 1) = cid;
                    }
                }
            }
            //remove currentP from seeds
            seeds.shed_row(0);
        }
    }
return clust;
}
