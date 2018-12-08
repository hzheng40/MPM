//
// Created by billyzheng on 10/24/18.
//

#ifndef MPM_CONSTANTS_H
#define MPM_CONSTANTS_H
#include <Eigen/Dense>
using namespace Eigen;

const double BSPLINE_EPSILON = 1e-4;
const double MIN_FLOAT = 1e-4;
const double MAX_FLOAT = 1e7;
const int BSPLINE_RADIUS = 2;

static const double
        DENSITY = 1e4,
        PARTICLE_DIAM = 500,
//        PT_MASS = 0.000001,
        PT_MASS = 1,
        CRIT_COMPRESS = 2.5e-2,
        CRIT_STRETCH = 7.5e-3,
        HARDENING = 10.0,
        YOUNGS_MODULUS = 50,
        POISSONS_RATIO = 0.3,
//        YOUNGS_MODULUS = 1.4e5,
//        POISSONS_RATIO = 0.2,
        STICKY = .4,
        BOUNCY = .9,
        GRAVITY = 9.8,
//        GRAVITY = 0,
//        GRAVITY = 9800,
        TIMESTEP = 5e-4,
        MAX_ITER = 20000;

static const double
        LAMBDA = YOUNGS_MODULUS*POISSONS_RATIO/((1+POISSONS_RATIO)*(1-2*POISSONS_RATIO)),
        MU = YOUNGS_MODULUS/(2+2*POISSONS_RATIO);

#endif //MPM_CONSTANTS_H
