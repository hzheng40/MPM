//
// Created by billyzheng on 10/24/18.
//

#ifndef MPM_CONSTANTS_H
#define MPM_CONSTANTS_H
#include <Eigen/Dense>
using namespace Eigen;

const float BSPLINE_EPSILON = 1e-4;
const int BSPLINE_RADIUS = 2;
static const int
        MAX_TIMESTEP = 1;
static const float
        PT_MASS = 0.1,
        CFL = .04,					//Adaptive timestep adjustment
//        MAX_TIMESTEP = 5e-4,		//Upper timestep limit?
        CRIT_COMPRESS = 1-1.9e-2,
        CRIT_STRETCH = 1+7.5e-3,
        HARDENING = 5.0,
        DENSITY = 100,
        YOUNGS_MODULUS = 1.5e5,
        POISSONS_RATIO = .2,
        STICKY = .9,				//Collision stickiness (lower = stickier)
        GRAVITY = -9.8,
        TIMESTEP = 0.0001;

//Actual timestep is adaptive, based on grid resolution and max velocity
//extern float TIMESTEP;

static const float
        LAMBDA = YOUNGS_MODULUS*POISSONS_RATIO/((1+POISSONS_RATIO)*(1-2*POISSONS_RATIO)),
        MU = YOUNGS_MODULUS/(2+2*POISSONS_RATIO);

#endif //MPM_CONSTANTS_H
