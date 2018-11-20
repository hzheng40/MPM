//
// Created by billyzheng on 10/24/18.
//

#include "Particle.h"
#include <iostream>
using namespace Eigen;
Particle::Particle() {}
Particle::Particle(Vector3f pos, Vector3f vel, float mass, float lambda, float mu, float timestep) {
    // load particle properties
    position = pos;
    velocity = vel;
    velocity_gradient = Matrix3f::Zero();
    this->mass = mass;
    this->lambda = lambda;
    this->mu = mu;
    this->timestep = timestep;
    def_elastic = Matrix3f::Identity();
    def_plastic = Matrix3f::Identity();

}
Particle::~Particle(){}
void Particle::updatePos() {
//    if (velocity(0)>1000.0) velocity(0) = 1000.0;
//    if (velocity(0)<-1000.0) velocity(0) = -1000.0;
//    if (velocity(1)>1000.0) velocity(1) = 1000.0;
//    if (velocity(1)<-1000.0) velocity(1) = -1000.0;
//    if (velocity(2)>1000.0) velocity(2) = 1000.0;
//    if (velocity(2)<-1000.0) velocity(2) = -1000.0;
    Vector3f change = timestep*velocity;
    if (isnan(change(0)) || isnan(change(1)) || isnan(change(2))) {
        return;
    }
    position += change;
}
void Particle::updateGradient() {
    velocity_gradient *= timestep;
    velocity_gradient += Matrix3f::Identity();
    def_elastic = velocity_gradient*def_elastic;
}
void Particle::applyPlasticity() {
    Matrix3f f_all = def_elastic * def_plastic;
    JacobiSVD<Matrix3f> svd(def_elastic, ComputeFullV | ComputeFullU);
    Vector3f val_singular = svd.singularValues();
    Matrix3f u = svd.matrixU();
    Matrix3f v = svd.matrixV();
    // clamp singular value to elastic
    for (int i = 0; i<3; i++) {
        float sig = val_singular(i);
        if (sig < CRIT_COMPRESS) val_singular(i) = CRIT_COMPRESS;
        else if (sig > CRIT_STRETCH) val_singular(i) = CRIT_STRETCH;
    }
    Matrix3f v_copy(v), u_copy(u);
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            v_copy(i,j) /= val_singular(i);
        }
    }
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            u_copy(i,j) *= val_singular(i);
        }
    }
    def_plastic = v_copy * u.transpose() * f_all;
    def_elastic = u_copy * v.transpose();
}

//const Matrix3f Particle::energyDerivative() {
//    float harden = exp(HARDENING*(1-def_plastic.determinant()));
//    float Je = svd_e(0)*svd_e(1)*svd_e(2);
//    Matrix3f corot = 2*mu*(def_elastic-svd_w*svd_v.transpose())*def_elastic.transpose();
//    corot += Matrix3f::Identity() * (lambda*Je*(Je-1));
//    return volume * harden * corot;
//}

// neo hookean
const Matrix3f Particle::energyDerivative() {
    JacobiSVD<Matrix3f> svd(def_elastic, ComputeFullV | ComputeFullU);
    Vector3f val_singular = svd.singularValues();
    Matrix3f u = svd.matrixU();
    Matrix3f v = svd.matrixV();
    // regular to polar svd
    if (u.determinant() == -1.0) {
        u(0, 2) = -u(0, 2); u(1, 2) = -u(1, 2); u(2, 2) = -u(2, 2);
        val_singular(2) = -val_singular(2);
    }
    if (v.determinant() == -1.0) {
        v(0, 2) = -v(0, 2); v(1, 2) = -v(1, 2); v(2, 2) = -v(2, 2);
        val_singular(2) = -val_singular(2);
    }
    // putting svd back together
    Matrix3f F = u*val_singular.asDiagonal()*v.transpose();
    float J = F.determinant();
//    float energy_density = (mu/2)*((F*F.transpose()).trace()-3) - mu*log(J) + (lambda/2)*pow((log(J)),2);
    Matrix3f P = mu*(F-F.transpose()) + lambda*log(J)*(F.inverse().transpose());
    Matrix3f energy = volume*P*def_elastic.transpose();
//    cout << "energy: \n" << energy << endl;
    return energy;
}
