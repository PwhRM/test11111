//
//  inImDynM.cpp
//  DSets-DBSCAN
//
//orignal paper:
// S. Rota Bulo, and I. M. Bomze. Infection and immunization: a new class of evolutionarygame dynamics.
//
//  Created by Rui Ma on 02.07.18.
//  Copyright © 2018 Rui Ma. All rights reserved.
//
//perform histogram equalization on input image
//inputs:
//1) mat A: a pairwise n * n similarity matrix (with zero diagonal)
//2) mat a: n * 1 arma matrix, n-dimensional standard simplex
//3) double toll: precision required from the dynamical system
//4) int maxIters: maximum number of iterations
//output:
//1) mat x: n * 1 arma matrix, the population vector at convergence
///////////////////////////////////////////////////////

#include "inImDynM.hpp"
#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

int selectPureStrategy(mat x, mat r);

mat inImDynM(mat A, mat x, double toll, int maxIters){
    int iters = 1;
    double NashError = 2 * toll;
    
    mat g = A * x;
    while (iters < maxIters && NashError >= toll) {
        mat xt = x.t() * g;
        mat r = g - xt(0, 0);
        //p-norm
        NashError = norm(arma::min(x, -r), 2);
        /*
        cout << iters << endl;
        if (iters == 277) {
            cout << iters << endl;
            vec testnash = arma::min(x, -r);
            cout << NashError << endl;
            cout << norm(testnash, 2) << endl;
        }
         */
        int i = selectPureStrategy(x, r);
        double den = A(i, i) - g(i, 0) - r(i, 0);
        bool do_remove = false;
        double mu;
        double optDelta;
        if (r(i, 0) >= 0) {
            mu = 1;
            if (den < 0) {
                mu = std::min(mu, - r(i) / den);
                if (mu < 0) {
                    mu = 0;
                }
            }
        } else {
            do_remove = true;
            mu = x(i, 0) / (x(i, 0) - 1);
            if (den < 0) {
                optDelta = - r(i, 0) / den;
                if (optDelta >= mu) {
                    mu = optDelta;
                    do_remove = false;
                }
            }
        }
        mat tmp = - x;
        tmp(i, 0) = tmp(i, 0) + 1;
        x = x + mu * tmp;
        if (do_remove) {
            x(i, 0) = 0;
        }
        x = abs(x) / arma::sum(arma::sum(abs(x)));
        g = mu * (A.col(i) - g) + g;
        iters = iters + 1;
    }
    cout << "number of iters in inImDynM: " << iters << endl;
    return x;
}

int selectPureStrategy(mat x, mat r){
    /*
    double maxVal = 0, minVal = 0;
    int i, maxIdx = -1, minIdx = -1;
    for (int j = 0; j < x.n_rows; j++) {
        if (maxVal < r(j, 0)) {
            maxVal = r(j, 0);
            maxIdx = j;
        } else if (minVal > r(j, 0) && x(j, 0) > 0) {
            minVal = r(j, 0);
            minIdx = j;
        }
    }
    i = maxVal >= -minVal ? maxIdx : minIdx;
    */
    uword i, j;
    //index
    vec index = zeros(x.n_elem);
    for (uword k = 0; k < x.n_elem; k++) {
        index(k) = double(k);
    }
    //mask_index
    vec masked_index = index(find(x > 0));
    //r_sub: r(x>0)
    //r is a n*1 mat so that can be used as vec
    vec r_sub = r(find(x > 0));
    i = x.index_max();
    j = r_sub.index_min();
    j = uword(masked_index(j));
    if (r(i) < - r(j)) {
        i = j;
    }
    return int(i);
}
