//
// Created by Hongrui Zheng on 10/23/18.
//

#include "MPM.h"

namespace MPM {
    Grid::Grid(int size, float resolution) {
        grid_size = size;
        grid_resolution = resolution;
        grid_ind = Eigen::MatrixXi::Zero(size, size);
        grid_mass = Eigen::MatrixXf::Zero(size, size);
        grid_vel = Eigen::MatrixXf::Zero(size, size);
    }
    Grid::~Grid() {
        std::cout << "Destroying Grid object." << std::endl;
    }

    void Grid::updateGrid(Particle& particle) {

    }
}