//
// Created by billyzheng on 10/24/18.
//

#include "Particle.h"
using namespace Eigen;
Particle::Particle() {}
Particle::Particle(const Vector2f &pos, const Vector2f &vel, float mass, float lambda, float mu, float timestep) {
    // load particle properties
    position = pos;
    velocity = vel;
    this->mass = mass;
    this->lambda = lambda;
    this->mu = mu;
    this->timestep = timestep;

}
Particle::~Particle(){}
void Particle::updatePos() {
    position += timestep*velocity;
}
void Particle::updateGradient() {
    velocity_gradient *= timestep;
    DiagonalMatrix<float, 2> m(1,1);
    velocity_gradient += m;
    def_elastic = velocity_gradient*def_elastic;
}
void Particle::applyPlasticity() {
    Matrix2f f_all = def_elastic * def_plastic;
    JacobiSVD<Matrix2f> svd(def_elastic, ComputeFullV | ComputeFullU);
//    def_elastic.svd(&svd_w, &svd_e, &svd_v);
}

const Matrix2f Particle::energyDerivative() {
//    float harden = exp();
}
