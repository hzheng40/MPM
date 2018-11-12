//
// Created by billyzheng on 11/11/18.
//

#ifndef PARTIO_EXAMPLE_POISSON_H
#define PARTIO_EXAMPLE_POISSON_H

#include <Eigen/Dense>
#include "Constants.h"
#include "Particle.h"
using namespace Eigen;
using namespace std;

class Poisson {
public:
    MatrixXd grid;
    int grid_width, grid_height;
    float cell_size;
    float radius;
    int k;

    Poisson(int grid_width, int grid_height, float radius, int k);
    virtual ~Poisson();
    void initGrid();
    void sample();
};
#endif //PARTIO_EXAMPLE_POISSON_H
