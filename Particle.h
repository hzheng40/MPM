//
// Created by billyzheng on 10/24/18.
//

#ifndef MPM_PARTICLE_H
#define MPM_PARTICLE_H
#include <Eigen/Dense>
#include <Eigen/SVD>
#include "Constants.h"
using namespace Eigen;
using namespace std;
class Particle {
public:
    double volume, mass, density;
    Vector3d position, velocity;
    Matrix3d velocity_gradient;
    double lambda, mu, timestep;
    Matrix3d def_elastic, def_plastic, full_def;
    Matrix3d svd_u, svd_v;
    Vector3d svd_s;
    Vector3d grid_position;
    // 64 in 3d (4x4x4)
    Vector3d weight_gradient[64];
    double weights[64];
    Particle();
    Particle(Vector3d pos, Vector3d vel, double mass, double lambda, double mu, double timestep);
    virtual ~Particle();
    // check weight sum
    void checksum();
    void checkgradsum();
    // update position based on velocity
    void updatePos();
    // update deformation gradient
    void updateGradient();
    void applyPlasticity();
    Matrix3d stress();

};

#endif //MPM_PARTICLE_H