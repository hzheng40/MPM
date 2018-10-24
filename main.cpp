#include <iostream>
#include "MPM.h"
int main() {
    std::cout << "Starting MPM big loop" << std::endl;
    int num_iter = 5000;
    MPM::Grid grid = MPM::Grid(500, 0.5);
    MPM::ParticleList particles = MPM::ParticleList();
    for (int i = 0; i < num_iter; ++i) {

    }
    return 0;
}