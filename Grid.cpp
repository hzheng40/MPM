//
// Created by billyzheng on 10/24/18.
//

#include "Grid.h"

Grid::Grid(Vector2f pos, Vector2f dims, Vector2f cells, Object* object) {
    obj = object;
    origin = pos;
    cellsize = dims/cells;
    size = cells+1;
    nodes_length = size[0]*size[1];
    nodes = new GridNode[nodes_length];
    node_area = cellsize[0]*cellsize[1];
}

Grid::~Grid() {
    delete[] nodes;
}

void Grid::initializeMass() {
    memset(nodes, 0, sizeof(GridNode)*nodes_length);
    for (int i = 0; i<obj->size)
}

void Grid::initilaizeVelocities() {

}

void Grid::calculateVolumes() const {

}

void Grid::explicitVelocities(const Vector2f &gravity) {
    for (int i=0; i<obj->size; i++) {
        Partile& p = obj->particles[i];
        Matrix2f energy = p.energyDerivative();
        int ox = p.grid_position[0], oy = p.grid_position[1];
        for (int idx=0, y=oy-1, y_end=y+3;  y<=y_end; y++) {
            for (int x=ox-1, x_end=x+3; x<=x_end; x++, idx++) {
                float w = p.weights[idx];
                if (w > BSPLINE_EPSILON) {
                    int n = (int) (y*size[0]+x);
                    nodes[n].force += energy*p.weight_gradient[idx];
                }
            }
        }
    }

    for (int i=0; i<nodes_length; i++) {
        GridNode &node = nodes[i];
        if (node.on) {
            node.force = node.velocity + TIMESTEP*(gravity - node.force/node.mass);
        }
    }
    collisionGrid();
}

void Grid::updateVelocities() const {
    for (int i=0; i<obj->size; i++) {
        Particle& p = obj->particles[i];
        Vector2f pic, flip = p.velocity;
        Matrix2f& grad = p.velocity_gradient;
        grad.setData(0.0);
        int ox = p.grid_position[0], oy = p.grid_position[1];
        for (int idx=0, y=oy-1, y_end=y+3; y<=y_end; y++) {
            for (int x=ox-1, x_end=x+3; x<=x_end; x++, idx++) {
                float w = p.weights[idx];
                if (w > BSPLINE_EPSILON) {
                    GridNode &node = nodes[(int) (y*size[0]+x)];
                    pic += node.force*w;
                    flip += (node.force - node.velocity)*w;
                    grad += node.force.outer_product(p.weight_gradient[idx]);
                }
            }
        }
    }
}