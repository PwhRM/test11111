//
//  dist_mat_maker.cpp
//  DSets-DBSCAN
//
//  Created by Rui Ma on 06.07.18.
//  Copyright © 2018 Rui Ma. All rights reserved.
//
//generate Euclidean distacne matrix
//input:
//1) mat data: arma matrix n * 2: col x, col y
//output:
//1) mat dist_mat: arma matrix: Euclidean distacne matrix

#include <armadillo>
#include <math.h>
#include "dist_mat_maker.hpp"

using namespace arma;
using namespace std;

mat dist_mat_maker(mat data){
    uword n_data = data.n_rows;
    mat dist_mat = zeros(n_data, n_data);
    
    for (uword i = 0; i < n_data; i++) {
        for (uword j = i; j < n_data; j++) {
            dist_mat(i, j) = sqrt(pow(data(i, 0) - data(j, 0), 2) + pow(data(i, 1) - data(j, 1), 2));
        }
    }
    
    for (uword i = 0; i < n_data; i++) {
        for (uword j = i; j < n_data; j++) {
            dist_mat(j, i) = dist_mat(i, j);
        }
    }
    
    return dist_mat;
}
