//
// Created by billyzheng on 10/24/18.
//

#ifndef MPM_GRID_H
#define MPM_GRID_H

#include <Eigen/Dense>
#include "Constants.h"
#include "Particle.h"
using namespace Eigen;
using namespace std;
typedef struct GridNode {
    // explicit
    float mass;
    bool on;
    Vector2f velocity, force;

} GridNode;

class Grid {
public:
    Vector2f origin, size, cellsize;
    array<Particle, TEST_SIZE> object;
    int nodes_length;
    float node_area, TIMESTEP;
    GridNode* nodes; //start of grid nodes

    Grid(Vector2f pos, Vector2f dims, Vector2f cells, array<Particle, TEST_SIZE> object);
    virtual ~Grid();
    // particles to grid
    void initializeMass();
    void initilaizeVelocities();
    void calculateVolumes() const;
    void calculateVelocities(const Vector2f &gravity);
    void updateVelocities();
    // collision
    void collisionGrid();
    void collisionParticles();
    // interpolation with cubic bspline
    static float bspline(float x) {
        x = fabs(x);
        float w;
        if (x < 1) {
            w = x*x*(x/2-1)+2/3.0;
        } else if (x < 2) {
            w = x*(x*(-x/6+1)-2)+4/3.0;
        } else {
            return 0;
        }
        if (w < BSPLINE_EPSILON) {
            return 0;
        }
        return w;
    }
    static float bsplinePrime(float x) {
        float abs_x = fabs(x), w;
        if (abs_x < 1) {
            return 1.5*x*abs_x-2*x;
        } else if (x < 2) {
            return -x*abs_x/2+2*x-2*x/abs_x;
        } else {
            return 0;
        }
    }
};
#endif //MPM_GRID_H
