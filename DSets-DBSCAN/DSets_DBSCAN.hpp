//
//  DSets_DBSCAN.hpp
//  DSets-DBSCAN
//
//  Created by Rui Ma on 29.06.18.
//  Copyright © 2018 Rui Ma. All rights reserved.
//

#ifndef DSets_DBSCAN_hpp
#define DSets_DBSCAN_hpp

#include <stdio.h>
#include <iostream>
#include <armadillo>
#include "hist_eq.hpp"
#include "inImDynM.hpp"

using namespace arma;

mat DSets_DBSCAN(mat dist_mat, double supportThreshold = 1e-4, double precision = 1e-8, int maxIters = 1000, int minPts = 4, int nbins = 50);
#endif /* DSets_DBSCAN_hpp */
