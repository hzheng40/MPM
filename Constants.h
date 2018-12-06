//
// Created by billyzheng on 10/24/18.
//

#ifndef MPM_CONSTANTS_H
#define MPM_CONSTANTS_H
#include <Eigen/Dense>
using namespace Eigen;

const float BSPLINE_EPSILON = 1e-7;
const float MIN_FLOAT = 1e-7;
const float MAX_FLOAT = 1e7;
const int BSPLINE_RADIUS = 3;

static const float
        DENSITY = 1,
        PARTICLE_DIAM = 500,
//        PT_MASS = 0.000001,
        PT_MASS = 1,
        CRIT_COMPRESS = 2.5e-2,
        CRIT_STRETCH = 7.5e-3,
        HARDENING = 10.0,
//        YOUNGS_MODULUS = 1.4e5,
//        POISSONS_RATIO = 0.2,
        YOUNGS_MODULUS = 1.4e5,
        POISSONS_RATIO = 0.2,
        STICKY = .4,
        BOUNCY = .9,
//        GRAVITY = 9.8,
        GRAVITY = 0,
//        GRAVITY = 980000,
        TIMESTEP = 1e-4,
        MAX_ITER = 300;

static const float
        LAMBDA = YOUNGS_MODULUS*POISSONS_RATIO/((1+POISSONS_RATIO)*(1-2*POISSONS_RATIO)),
        MU = YOUNGS_MODULUS/(2+2*POISSONS_RATIO);

#endif //MPM_CONSTANTS_H
