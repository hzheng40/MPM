//
// Created by billyzheng on 10/24/18.
//

#ifndef MPM_PARTICLE_H
#define MPM_PARTICLE_H
#include <Eigen/Dense>

namespace MPM {
    class Particle {
    public:
        float volume, mass, density;
        Eigen::Vector2f position, velocity;
        Eigen::Matrix2f velocity_gradient;
        // lame params
        float lambda, mu;
        // deformation gradient
        Eigen::Matrix2f def_elastic, def_plastic;
        // cached svd's for elastic deformation gradient
        Eigen::Matrix2f svd_w, svd_v;
        Eigen::Vector2f svd_e;
        // cached polar decomp
        Eigen::Matrix2f polar_r, polar_s;
        // grid interpolation weights
        Eigen::Vector2f grid_position;
        Eigen::Vector2f weight_gradient[16];
        float weights[16];

        Particle();
        Particle(const Eigen::Vector2f& pos, const Eigen::Vector2f& vel, float mass, float lambda, float mu);
        virtual ~Particle();

        // update position based on velocity
        void updatePos();
        // update deformation gradient
        void updateGradient();
        void applyPlasticity();
        // compute stress tensor
        const Eigen::Matrix2f energyDerivative();
        // compute stress force delta, implicit velocity update
        const Eigen::Vector2f deltaForce(const Eigen::Vector2f& u, const Eigen::Vector2f& weight_grad);
    };
}

#endif //MPM_PARTICLE_H