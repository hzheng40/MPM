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
//    this->volume = this->mass/DENSITY;
    this->mu = mu;
    this->timestep = timestep;
    def_elastic = Matrix3f::Identity();
    def_plastic = Matrix3f::Identity();

}
Particle::~Particle(){}

void Particle::checksum() {
    float sum;
    for (int i=0; i<64; i++) {
        sum += weights[i];
    }
    cout << "weight sum: " << sum << endl;
}

void Particle::checkgradsum() {
    float sum;
    for (int i=0; i<64; i++) {
        sum += weight_gradient[i](0) + weight_gradient[i](1) + weight_gradient[i](2);
    }
    cout << "weight grad sum: " << sum << endl;
}

void Particle::updatePos() {
//    if (velocity(0)>MAX_FLOAT) velocity(0) = MAX_FLOAT;
//    if (fabs(velocity(0))<MIN_FLOAT) velocity(0) = velocity(0)/fabs(velocity(0)) * MIN_FLOAT;
//    if (velocity(1)>MAX_FLOAT) velocity(1) = MAX_FLOAT;
//    if (fabs(velocity(1))<MIN_FLOAT) velocity(1) = velocity(1)/fabs(velocity(1)) * MIN_FLOAT;
//    if (velocity(2)>MAX_FLOAT) velocity(2) = MAX_FLOAT;
//    if (fabs(velocity(2))<MIN_FLOAT) velocity(2) = velocity(2)/fabs(velocity(2)) * MIN_FLOAT;

    Vector3f change = timestep*velocity;
//    if (isnan(change(0)) || isnan(change(1)) || isnan(change(2))) {
//        return;
//    }
    position += change;
}
void Particle::updateGradient() {
//    velocity_gradient *= timestep;
//    Matrix3f gradient = velocity_gradient + Matrix3f::Identity();
    def_elastic = def_elastic + timestep*velocity_gradient*def_elastic;
//    cout << "def_elastic: \n" << def_elastic << endl;
//    full_def = full_def + timestep*velocity_gradient*full_def;
}
void Particle::applyPlasticity() {
    def_elastic = svd_u * svd_s.asDiagonal() * svd_v.transpose();
    def_plastic = svd_v * svd_s.asDiagonal().inverse() * svd_u.transpose() * full_def;
}

//const Matrix3f Particle::energyDerivative() {
//    float harden = exp(HARDENING*(1-def_plastic.determinant()));
//    float Je = svd_s(0)*svd_s(1)*svd_s(2);
//    Matrix3f corot = 2*mu*(def_elastic-svd_u*svd_v.transpose())*def_elastic.transpose();
//    corot += Matrix3f::Identity() * (lambda*Je*(Je-1));
//    return this->volume * harden * corot;
//}

// neo hookean
const Matrix3f Particle::energyDerivative() {
    JacobiSVD<Matrix3f> svd(def_elastic, ComputeFullV | ComputeFullU);
    Vector3f val_singular = svd.singularValues();
    Matrix3f u = svd.matrixU();
    Matrix3f v = svd.matrixV();
//    cout << val_singular << endl;
//    cout << u << endl;
//    cout << v << endl;
    // regular to polar svd
//    if (u.determinant() == -1.0) {
    if (u.determinant() < 0.0) {
        u(0, 2) = -u(0, 2); u(1, 2) = -u(1, 2); u(2, 2) = -u(2, 2);
        val_singular(2) = -val_singular(2);
    }
//    if (v.determinant() == -1.0) {
    if (v.determinant() < 0.0) {
        v(0, 2) = -v(0, 2); v(1, 2) = -v(1, 2); v(2, 2) = -v(2, 2);
        val_singular(2) = -val_singular(2);
    }

    if (val_singular(0) < val_singular(1)) {
        float temp = val_singular(0);
        val_singular(0) = val_singular(1);
        val_singular(1) = temp;
    }
    for (int i = 0; i<3; i++) {
//        float sig = val_singular(i);
        val_singular(i) = max(1-CRIT_COMPRESS, min(val_singular(1), 1+CRIT_STRETCH));
//        if (sig < CRIT_COMPRESS) val_singular(i) = CRIT_COMPRESS;
//        else if (sig > CRIT_STRETCH) val_singular(i) = CRIT_STRETCH;
    }
    svd_u = Matrix3f(u);
    svd_v = Matrix3f(v);
    svd_s = Vector3f(val_singular);
    // putting svd back together

//    def_elastic = u*val_singular.asDiagonal()*v.transpose();

//    float J = def_elastic.determinant();
//    float energy_density = (mu/2)*((F*F.transpose()).trace()-3) - mu*log(J) + (lambda/2)*pow((log(J)),2);
    // neo-Hookean
//    Matrix3f F = u*val_singular.asDiagonal()*v.transpose();
//    float J = F.determinant();
//    Matrix3f P = mu*(F-F.inverse().transpose()) + lambda*log(J)*(F.inverse().transpose());
    // corot
    Matrix3f F = u*v.transpose()*v*val_singular.asDiagonal()*v.transpose();
    float J = F.determinant();
    Matrix3f R = u*v.transpose();
    Matrix3f P = 2*mu*(F-R) + lambda*(J-1)*J*(F.inverse().transpose());
//    Matrix3f P = mu*(def_elastic-def_elastic.transpose()) + lambda*log(J)*(def_elastic.inverse().transpose());
    Matrix3f energy = this->volume*P*F.transpose();
//    Matrix3f energy = 10000000*P*def_elastic.transpose();
//    cout << "energy: \n" << energy << endl;
//    cout << "def_elastic: \n" << def_elastic << endl;
//    cout << "volume: " << this->volume << endl;
    return energy;
//    return Matrix3f::Identity();
}

