//
// Created by billyzheng on 10/24/18.
//

#include "Particle.h"
using namespace Eigen;
Particle::Particle() {}
Particle::Particle(Vector3f pos, Vector3f vel, float mass, float lambda, float mu, float timestep) {
    // load particle properties
    position = pos;
    velocity = vel;
    this->mass = mass;
    this->lambda = lambda;
    this->mu = mu;
    this->timestep = timestep;
    def_elastic = Matrix3f::Identity();
    def_plastic = Matrix3f::Identity();

}
Particle::~Particle(){}
void Particle::updatePos() {
    position += timestep*velocity;
}
void Particle::updateGradient() {
    velocity_gradient *= timestep;
    velocity_gradient += Matrix3f::Identity();
    def_elastic = velocity_gradient*def_elastic;
}
void Particle::applyPlasticity() {
    Matrix3f f_all = def_elastic * def_plastic;
    JacobiSVD<Matrix3f> svd(def_elastic, ComputeFullV | ComputeFullU);
}

const Matrix3f Particle::energyDerivative() {
    float harden = exp(HARDENING*(1-def_plastic.determinant()));
    float Je = svd_e(0)*svd_e(1)*svd_e(2);
    Matrix3f corot = 2*mu*(def_elastic-svd_w*svd_v.transpose())*def_elastic.transpose();
    corot += Matrix3f::Identity() * (lambda*Je*(Je-1));
    return volume * harden * corot;
}
