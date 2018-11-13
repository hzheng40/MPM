//
// Created by billyzheng on 11/11/18.
//

#ifndef PARTIO_EXAMPLE_POISSON_H
#define PARTIO_EXAMPLE_POISSON_H

#define EIGEN_USE_NEW_STDVECTOR

#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <vector>
#include <math.h>
#include "Constants.h"
#include "Particle.h"
using namespace Eigen;
using namespace std;

class Poisson {
public:
    vector<Vector3f, aligned_allocator<Vector3f>> grid;
    // indexing: i*w + j + k*w*h
    // i: y, j: x, k: z
    int grid_width, grid_height, grid_depth;
    float cell_size;
    float radius;
    int k;
    vector<int> active;

    Poisson(int grid_width, int grid_height, int grid_depth, float radius, int k);
    virtual ~Poisson();
    void initGrid();
    void sample();
    void writePartio(const string& particle_file);
    void writePartioByFrame(const string& particle_file_prefix, int seq_num);
};
#endif //PARTIO_EXAMPLE_POISSON_H
