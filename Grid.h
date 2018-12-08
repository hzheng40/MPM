//
// Created by billyzheng on 10/24/18.
//

#ifndef MPM_GRID_H
#define MPM_GRID_H

#include <Eigen/Dense>
#include <vector>
#include "Constants.h"
#include "Particle.h"
using namespace Eigen;
using namespace std;
typedef struct GridNode {
    // explicit
    double mass;
    bool on;
    Vector3d velocity, grid_forces;
} GridNode;

class Grid {
public:
    Vector3d origin, size, cellsize;
    vector<Particle> object;
    int nodes_length;
    double node_volume, timestep;
    GridNode* nodes; //start of grid nodes

    Grid(Vector3d pos, Vector3d dims, Vector3d cell_num, vector<Particle> object);
    virtual ~Grid();
    // particles to grid
    void initializeMass();
    void initializeVel();
    void calculateVolumes();
    void p2g_vel(const Vector3d &gravity);
    void g2p_vel();
    // collision
    void collisionGrid();
    void collisionParticles();
    // interpolation with cubic bspline
    static double bspline(double x) {
        x = std::abs(x);
        double w;
        if (x < 1) {
            w = x*x*(x/2-1)+2/3.0;
        } else if (x < 2) {
            w = x*(x*(-x/6+1)-2)+4/3.0;
        } else {
            return 0.0;
        }
        if (w < BSPLINE_EPSILON) {
            return 0.0;
        }
        return w;
    }
    static double bsplinePrime(double x) {
        double abs_x = std::abs(x), w;
        if (abs_x < 1) {
            return 1.5*abs_x*abs_x-2*x;
        } else if (x < 2) {
            return -x*abs_x/2+2*x-2*x/abs_x;
        } else {
            return 0.0;
        }
    }
};
#endif //MPM_GRID_H
