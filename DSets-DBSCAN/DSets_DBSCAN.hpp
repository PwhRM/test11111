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

mat DSets_DBSCAN(mat dist_mat, vec x0, double supportThreshold, double precision, int maxIters, int minPts);
#endif /* DSets_DBSCAN_hpp */
