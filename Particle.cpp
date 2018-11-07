//
// Created by billyzheng on 10/24/18.
//

#include "Particle.h"
using namespace Eigen;
Particle::Particle() {}
Particle::Particle(const Vector2f &pos, const Vector2f &vel, float mass, float lambda, float mu) {
    // load particle properties
    position = pos;
    velocity = vel;
    this->mass = mass;
    this->lambda = lambda;
    this->mu = mu;

}
Particle::~Particle(){}
void Particle::updatePos() {
    position += TIMESTEP*velocity;
}
void Particle::updateGradient() {
    velocity_gradient *= TIMESTEP;
    DiagonalMatrix<float, 2> m(1,1);
    velocity_gradient += m;
    def_elastic = velocity_gradient*def_elastic;
}
void Particle::applyPlasticity() {

}
