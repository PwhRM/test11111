//
//  inImDynM.hpp
//  DSets-DBSCAN
//
//  Created by Rui Ma on 02.07.18.
//  Copyright © 2018 Rui Ma. All rights reserved.
//

#ifndef inImDynM_hpp
#define inImDynM_hpp

#include <stdio.h>
#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

mat inImDynM(mat A, mat x, double toll, int maxIters);
#endif /* inImDynM_hpp */
