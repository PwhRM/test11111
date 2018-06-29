//
//  hist_eq.cpp
//  DSets-DBSCAN
//
//  Created by Rui Ma on 29.06.18.
//  Copyright © 2018 Rui Ma. All rights reserved.
//
//perform histogram equalization on input image
//inputs:
//1) mat img: arma mat matrix, input image
//2) int hist_size: number of bins in histogram
//output:
//1) mat new_img: arma mat matrix, output image
///////////////////////////////////////////////////////

#include "hist_eq.hpp"
#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

mat hist_eq(mat img, int hist_size){
    mat new_img;
    double max_value = img.max();
    //calculate cdf
    mat cdf(hist_size, 1);
    for (int i = 0; i < int(img.n_rows); i++){
        for(int j = 0; j < int(img.n_cols); j++){
            cdf(floor(img(i, j) / (max_value / hist_size)), 1)++;
        }
    }
    //apply transformation
    for (int i = 0; i < int(img.n_rows); i++){
        for(int j = 0; j < int(img.n_cols); j++){
            new_img(i, j) = cdf(img(i, j), 1) / (hist_size - 1);
        }
    }
    
    return new_img;
}
