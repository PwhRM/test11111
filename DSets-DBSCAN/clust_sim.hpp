//
//  clust_sim.hpp
//  DSets-DBSCAN
//
//  Created by Rui Ma on 29.06.18.
//  Copyright © 2018 Rui Ma. All rights reserved.
//

#ifndef clust_sim_hpp
#define clust_sim_hpp

#include <stdio.h>
#include "clust_sim.hpp"
#include <iostream>
#include <math.h>
#include <armadillo>

mat clust_sim(int n_clust, int n_point, int n_noise, double r, double height, double width);

#endif /* clust_sim_hpp */
