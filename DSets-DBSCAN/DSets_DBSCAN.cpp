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
//5) int maxIters: maximum number of iterations of the dynamic systems (default 1000)
//output:
//1) mat clust: arma mat matrix, clust id
///////////////////////////////////////////////////////
#include "DSets_DBSCAN.hpp"
#include <iostream>
#include <armadillo>
#include "hist_eq.hpp"
#include "inImDynM.hpp"
#include "ExpandCluster.hpp"

using namespace std;
using namespace arma;

vec ExpandCluster(mat dist_mat, uword pt, int cid, double Eps, int minPts, vec clust, vec d);

mat DSets_DBSCAN(mat dist_mat, double supportThreshold, double precision, int maxIters, int minPts, int nbins){
    //control of input arguments be done in main function
    //generate similarity matrix
    double sigma = 10 * mean(mean(dist_mat));
    mat A = exp(- dist_mat / sigma);
    //dignal should be zeros
    A.diag() = A.diag() - A.diag();
    //perform histogram equalization on similarity matrix
    A = hist_eq(A, nbins);
    //number of points in the data set
    int n_points = int(A.n_rows);
    //index of points in the orignal data set
    //d will be changing over the while loop
    //0 value will be occupied as visited points so the index starts from 1 instead
    vec d(A.n_rows);
    for (int i = 0; i < n_points; i++) {
        d(i) = i + 1;
    }
    //initial cluster id
    int cid = 0;
    //cluster results
    vec clust = zeros(n_points) - 1;
    //starting clustering
    //the visited data points will be assigned with 0 in the d vec after the while loop
    //first initialize count as n_points
    int count = n_points;
    while (any(d) && count > 0) {
        cid++;
        //frist remove visited data points
        vec d_tmp = d(find(d > 0));
        int n_points_tmp = int(d_tmp.n_elem);
        count = n_points_tmp;
        mat A_tmp(n_points_tmp, n_points_tmp),
        dist_mat_tmp(n_points_tmp, n_points_tmp);
        for (int i = 0; i < n_points_tmp; i++) {
            for (int j = 0; j < n_points_tmp; j++) {
                A_tmp(i, j) = A(d_tmp(i) - 1, d_tmp(j) - 1);
                dist_mat_tmp(i, j) = dist_mat(d_tmp(i) - 1, d_tmp(j) - 1);
            }
        }
        //get current x from current d
        //the size of d is changing, so size of x is changing too
        mat x(d_tmp.n_elem, 1);
        x = ones(d_tmp.n_elem, 1);
        x = x / arma::sum(arma::sum(x));
        //extract dominant from DSets method
        x = inImDynM(A_tmp, x, precision, maxIters);
        //note id_A will be vec instead of matrix
        //x_good_ind: index of good x in vector x
        uvec x_good_ind = find(x > supportThreshold);
        uvec id_A = x_good_ind;
        //assign clust id
        //id_clust - 1 is the actual index in clust
        vec id_clust = d_tmp.elem(x_good_ind);
        for (int i = 0; i < int(id_clust.n_elem); i++) {
            clust(id_clust(i) - 1) = cid;
        }
        //use DBSCAN to expand dominant
        double Eps = 0;
        for (int i = 0; i < int(id_A.n_elem); i++) {
            //dist_tmp should be the distance map in dist_mat, not dist_mat_tmp
            vec dist_tmp = dist_mat_tmp.col(id_A(i));
            dist_tmp = arma::sort(dist_tmp);
            Eps = max(Eps, dist_tmp(std::min(minPts, int(dist_tmp.n_elem)) - 1));
        }
        uvec visit_list = x_good_ind;
        for (int i = 0; i < int(x_good_ind.n_elem); i++) {
            uword pt = visit_list(i);
            clust = ExpandCluster(dist_mat_tmp, pt, cid, Eps, minPts, clust, d_tmp);
        }
        //get new xid
        vec d_clust = d(find(clust == cid));
        //assigan visited data points in d as 0
        for (int i = 0; i < d_clust.n_elem; i++) {
            d(d_clust(i) - 1) = 0;
        }
    }
    
    return clust;
}
