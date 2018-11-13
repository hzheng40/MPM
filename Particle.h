//
// Created by billyzheng on 10/24/18.
//

#ifndef MPM_PARTICLE_H
#define MPM_PARTICLE_H
#include <Eigen/Dense>
#include <Eigen/SVD>
#include "Constants.h"
using namespace Eigen;
class Particle {
public:
    float volume, mass, density;
    Vector2f position, velocity;
    Matrix2f velocity_gradient;
    float lambda, mu, timestep;
    Matrix2f def_elastic, def_plastic;
    Matrix2f svd_w, svd_v;
    Vector2f svd_e;
    Matrix2f polar_r, polar_s;
    Vector2f grid_position;
    // 64 in 3d (4x4x4)
    Vector3f weight_gradient[64];
    float weights[64];
    Particle();
    Particle(const Vector2f& pos, const Vector2f& vel, float mass, float lambda, float mu, float timestep);
    virtual ~Particle();
    // update position based on velocity
    void updatePos();
    // update deformation gradient
    void updateGradient();
    void applyPlasticity();
    // compute stress tensor
    const Matrix2f energyDerivative();
    // compute stress force delta, implicit velocity update
    const Vector2f deltaForce(const Vector2f& u, const Vector2f& weight_grad);
};

#endif //MPM_PARTICLE_H