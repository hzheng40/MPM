//
// Created by billyzheng on 10/24/18.
//

#include "Particle.h"
#include <iostream>
using namespace Eigen;
Particle::Particle() {}
Particle::Particle(Vector3d pos, Vector3d vel, double mass, double lambda, double mu, double timestep) {
    // load particle properties
    position = pos;
    velocity = vel;
    velocity_gradient = Matrix3d::Zero();
    this->mass = mass;
    this->lambda = lambda;
//    this->volume = this->mass/DENSITY;
    this->mu = mu;
    this->timestep = timestep;
    def_elastic = 1.8 * Matrix3d::Identity();
    def_plastic = 1.8 * Matrix3d::Identity();

}
Particle::~Particle(){}

void Particle::checksum() {
    double sum;
    for (int i=0; i<64; i++) {
        sum += weights[i];
    }
    cout << "weight sum: " << sum << endl;
}

void Particle::checkgradsum() {
    double sum;
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

    Vector3d change = timestep*velocity;
//    if (isnan(change(0)) || isnan(change(1)) || isnan(change(2))) {
//        velocity = Vector3d::Zero();
//        return;
//    }
    position += change;
}
void Particle::updateGradient() {
//    velocity_gradient *= timestep;
//    Matrix3f gradient = velocity_gradient + Matrix3f::Identity();
    Matrix3d change = timestep*velocity_gradient*def_elastic;
    if (isnan(change.sum())){
//        cout<<"change too small" << endl;
        return;
    }
    def_elastic = def_elastic + change;
//    cout << "defgrad: \n" << def_elastic << endl;
//    full_def = full_def + timestep*velocity_gradient*full_def;
}
void Particle::applyPlasticity() {
    JacobiSVD<Matrix3d> svd(def_elastic, ComputeFullV | ComputeFullU);
    Vector3d val_singular = svd.singularValues();
    Matrix3d u = svd.matrixU();
    Matrix3d v = svd.matrixV();
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
        double temp = val_singular(0);
        val_singular(0) = val_singular(1);
        val_singular(1) = temp;
    }
    for (int i = 0; i<3; i++) {
//        float sig = val_singular(i);
        val_singular(i) = max(1-CRIT_COMPRESS, min(val_singular(1), 1+CRIT_STRETCH));
//        if (sig < CRIT_COMPRESS) val_singular(i) = CRIT_COMPRESS;
//        else if (sig > CRIT_STRETCH) val_singular(i) = CRIT_STRETCH;
    }
    this->svd_u = Matrix3d(u);
    this->svd_v = Matrix3d(v);
    this->svd_s = Vector3d(val_singular);
    // putting svd back together

    this->def_elastic = u*val_singular.asDiagonal()*v.transpose();
    this->def_plastic = svd_v * svd_s.asDiagonal().inverse() * svd_u.transpose() * full_def;
}

//const Matrix3f Particle::stress() {
//    JacobiSVD<Matrix3f> svd(def_elastic, ComputeFullV | ComputeFullU);
//    Vector3f val_singular = svd.singularValues();
//    Matrix3f u = svd.matrixU();
//    Matrix3f v = svd.matrixV();
//    //    cout << val_singular << endl;
//    //    cout << u << endl;
//    //    cout << v << endl;
//    // regular to polar svd
//    //    if (u.determinant() == -1.0) {
//    if (u.determinant() < 0.0) {
//    u(0, 2) = -u(0, 2); u(1, 2) = -u(1, 2); u(2, 2) = -u(2, 2);
//    val_singular(2) = -val_singular(2);
//    }
//    //    if (v.determinant() == -1.0) {
//    if (v.determinant() < 0.0) {
//    v(0, 2) = -v(0, 2); v(1, 2) = -v(1, 2); v(2, 2) = -v(2, 2);
//    val_singular(2) = -val_singular(2);
//    }
//
//    if (val_singular(0) < val_singular(1)) {
//    float temp = val_singular(0);
//    val_singular(0) = val_singular(1);
//    val_singular(1) = temp;
//    }
//    for (int i = 0; i<3; i++) {
//    //        float sig = val_singular(i);
//    val_singular(i) = max(1-CRIT_COMPRESS, min(val_singular(1), 1+CRIT_STRETCH));
//    //        if (sig < CRIT_COMPRESS) val_singular(i) = CRIT_COMPRESS;
//    //        else if (sig > CRIT_STRETCH) val_singular(i) = CRIT_STRETCH;
//    }
//    svd_u = Matrix3f(u);
//    svd_v = Matrix3f(v);
//    svd_s = Vector3f(val_singular);
//    float harden = exp(HARDENING*(1-def_plastic.determinant()));
//    float Je = svd_s(0)*svd_s(1)*svd_s(2);
//    Matrix3f corot = 2*mu*(def_elastic-svd_u*svd_v.transpose())*def_elastic.transpose();
//    corot += Matrix3f::Identity() * (lambda*Je*(Je-1));
//    return this->volume * harden * corot;
//}

Matrix3d Particle::stress() {
    double harden = exp(HARDENING*(1-def_plastic.determinant()));

    // neo-Hookean
    double J = def_elastic.determinant();
    Matrix3d P = mu*(def_elastic - def_elastic.transpose().inverse()) + lambda*log(J)*def_elastic.transpose().inverse();
//    Matrix3d energy = this->volume*P*def_elastic.transpose();
//    Matrix3d stress = P*def_elastic.transpose();

    // corot
//    double J = def_elastic.determinant();
//    Matrix3d R = svd_u*svd_v.transpose();
//    Matrix3d P = 2*mu*(def_elastic-R) + lambda*(J-1)*J*(def_elastic.inverse().transpose());
//    Matrix3d energy = this->volume*harden*P;

//    cout << "stress: \n" << P << endl;
//    cout << "def_elastic: \n" << def_elastic << endl;
//    cout << "volume: " << this->volume << flush;
    return P;
//    return Matrix3f::Identity();
}

