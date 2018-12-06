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
#include "Partio.h"
#include <random>
#include <iostream>
using namespace Eigen;
using namespace std;

class Poisson {
public:
    vector<Vector3f, aligned_allocator<Vector3f>> grid;
    // indexing: i*w + j + k*w*h
    // i: y, j: x, k: z
    float grid_width, grid_height, grid_depth;
    int grid_length, w_ind, h_ind, d_ind;
    float cell_size;
    float radius;
    int k;
    vector<int> active;

    Poisson(int grid_width, int grid_height, int grid_depth, float radius, int k);
    virtual ~Poisson();
    void initGrid();
    void sample();
    vector<Vector3f, aligned_allocator<Vector3f>> toSphere();
    vector<Vector3f, aligned_allocator<Vector3f>> toCube();
    vector<Particle> toObject(vector<Vector3f, aligned_allocator<Vector3f>>);
    void writePartio(const string& particle_file);
    void writePartioByObejct(const string& particle_file, vector<Vector3f, aligned_allocator<Vector3f>> object_file);
    void writePartioByFrame(const string& particle_file_prefix, int seq_num);
};
#endif //PARTIO_EXAMPLE_POISSON_H
