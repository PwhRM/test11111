//
//  ExpandCluster.hpp
//  DSets-DBSCAN
//
//  Created by Rui Ma on 05.07.18.
//  Copyright © 2018 Rui Ma. All rights reserved.
//

#ifndef ExpandCluster_hpp
#define ExpandCluster_hpp

#include <stdio.h>
#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

vec ExpandCluster(mat dist_mat, uword pt, int cid, double Eps, int minPts, vec clust, vec d);
#endif /* ExpandCluster_hpp */
