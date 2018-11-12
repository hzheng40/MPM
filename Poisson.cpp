//
// Created by billyzheng on 11/11/18.
//

#include "Poisson.h"
using namespace Eigen;
Poisson::Poisson(int grid_width, int grid_height, float radius, int k) {
    //grid = MatrixXd(grid_height, grid_width);
    this->grid_width = grid_width;
    this->grid_height = grid_height;
    this->radius = radius;
    this->k = k;
    cell_size = radius/sqrt(2.0);
}

Poisson::~Poisson() {
}

void Poisson::initGrid() {
    grid = MatrixXd::Constant(grid_height, grid_width, -1.0);

}