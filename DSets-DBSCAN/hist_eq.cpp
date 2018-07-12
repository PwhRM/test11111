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
    mat new_img(img.n_rows, img.n_cols);
    double max_value = img.max();
    //normalize input image
    img = img / max_value;
    double bin_size = 1 / double(hist_size - 1);
    //calculate cdf
    mat cdf = zeros(hist_size , 1);
    mat pdf = zeros(hist_size, 1);
    for (int i = 0; i < int(img.n_rows); i++){
        for(int j = 0; j < int(img.n_cols); j++){
            pdf(floor(img(i, j) / bin_size), 0) =
            pdf(floor(img(i, j) / bin_size), 0) + 1;
        }
    }
    pdf = pdf / img.n_elem;
    for (uword i = 0; i < pdf.n_elem; i++) {
        for (uword j = 0; j <= i; j++) {
            cdf(i) = cdf(i) + pdf(j);
        }
    }
    //apply transformation
    //note that the value range for the cdf should be the same with the gray scale max of the image
    for (int i = 0; i < int(img.n_rows); i++){
        for(int j = 0; j < int(img.n_cols); j++){
            new_img(i, j) = cdf(floor(img(i, j) / bin_size), 0) * max_value;
        }
    }
    
    return new_img;
}
