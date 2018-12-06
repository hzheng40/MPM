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
    float volume, mass, density;
    Vector3f position, velocity;
    Matrix3f velocity_gradient;
    float lambda, mu, timestep;
    Matrix3f def_elastic, def_plastic, full_def;
    Matrix3f svd_u, svd_v;
    Vector3f svd_s;
    Vector3f grid_position;
    // 64 in 3d (4x4x4)
    Vector3f weight_gradient[64];
    float weights[64];
    Particle();
    Particle(Vector3f pos, Vector3f vel, float mass, float lambda, float mu, float timestep);
    virtual ~Particle();
    // check weight sum
    void checksum();
    void checkgradsum();
    // update position based on velocity
    void updatePos();
    // update deformation gradient
    void updateGradient();
    void applyPlasticity();
    const Matrix3f energyDerivative();

};

#endif //MPM_PARTICLE_H