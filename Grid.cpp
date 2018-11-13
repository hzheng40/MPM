//
// Created by billyzheng on 10/24/18.
//

#include "Grid.h"
Grid::Grid(Vector3f pos, Vector3f dims, Vector3f cell_num, vector<Particle> object) {
    origin = Vector3f();
    origin(0) = pos(0);
    origin(1) = pos(1);
    origin(2) = pos(2);
    cellsize = Vector3f();
    cellsize(0) = dims(0)/cell_num(0);
    cellsize(1) = dims(1)/cell_num(1);
    cellsize(2) = dims(2)/cell_num(2);
    size = cell_num + Vector3f::Ones();
    nodes_length = size(0)*size(1)*size(2);
    nodes = new GridNode[nodes_length];
    node_volume = cellsize(0)*cellsize(1)*cellsize(2);
    this->object = object;
    // timestep depends on grid resolution?
    timestep = 0.01;
}

Grid::~Grid() {
    delete[] nodes;
}

// initializing grid mass from original particles, only ran once at the beginning
void Grid::initializeMass() {
    memset(nodes, 0, sizeof(GridNode)*nodes_length);
    for (int i=0; i<object.size(); i++) {
        Particle& p = object[i];
        Vector3f particle_dist = p.position - origin;
        p.grid_position(0) = particle_dist(0)/cellsize(0);
        p.grid_position(1) = particle_dist(1)/cellsize(1);
        p.grid_position(2) = particle_dist(2)/cellsize(2);
        float px_grid = static_cast<int>(p.grid_position(0)),
              py_grid = static_cast<int>(p.grid_position(1)),
              pz_grid = static_cast<int>(p.grid_position(2));

        for (int idx = 0, x_offset=-1; x_offset <= 2; x_offset++) {

            for (int y_offset=-1; y_offset <= 2; y_offset++) {

                for (int z_offset=-1; z_offset <= 2; z_offset++, idx++) {

                    
                }
            }
        }




        float ox = p.grid_position(0), oy = p.grid_position(1);
        for (int idx=0, y=oy-1, y_end=y+3; y<=y_end; y++) {
            float y_pos = oy-y, wy = Grid::bspline(y_pos), dy = Grid::bsplinePrime(y_pos);
            for (int x=ox-1, x_end=x+3; x<=x_end; x++, idx++) {
                float x_pos = ox-x, wx = Grid::bspline(x_pos), dx = Grid::bsplinePrime(x_pos);
                float weight = wx*wy;
                p.weights[idx] = weight;
                p.weight_gradient[idx] = Vector2f(dx*wy, wx*dy);
                p.weight_gradient[idx](0) /= cellsize(0);
                p.weight_gradient[idx](1) /= cellsize(1);
                nodes[(int) (y*size(0)+x)].mass += weight*p.mass;
            }
        }
    }
}
// interpolate vel after mass, conserving momentum
void Grid::initilaizeVelocities() {
    for (int i=0; i<object.size(); i++) {
        Particle& p = object[i];
        float ox = p.grid_position(0), oy = p.grid_position(1);
        for (int idx=0, y=oy-1, y_end=y+3; y<=y_end; y++) {
            for (int x=ox-1, x_end=x+3; x<=x_end; x++, idx++) {
                float w = p.weights[idx];
                if (w > BSPLINE_EPSILON) {
                    int n = (int) (y*size(0)+x);
                    nodes[n].velocity += p.velocity * w * p.mass;
                    nodes[n].on = true;
                }
            }
        }
    }
    for (int i=0; i<nodes_length; i++) {
        GridNode &node = nodes[i];
        if (node.on) {
            node.velocity /= node.mass;
        }
    }
    collisionGrid();
}

void Grid::calculateVolumes() const {

}

void Grid::calculateVelocities(const Vector2f &gravity) {
    for (int i=0; i<object.size(); i++) {
        Particle& p = object[i];
        Matrix2f energy = p.energyDerivative();
        int ox = p.grid_position[0], oy = p.grid_position[1];
        for (int idx=0, y=oy-1, y_end=y+3;  y<=y_end; y++) {
            for (int x=ox-1, x_end=x+3; x<=x_end; x++, idx++) {
                float w = p.weights[idx];
                if (w > BSPLINE_EPSILON) {
                    int n = (int) (y*size(0)+x);
                    nodes[n].force += energy*p.weight_gradient[idx];
                }
            }
        }
    }

    for (int i=0; i<nodes_length; i++) {
        GridNode &node = nodes[i];
        if (node.on) {
            node.force = node.velocity + timestep*(gravity - node.force/node.mass);
        }
    }
    collisionGrid();
}

void Grid::updateVelocities() {
    for (int i=0; i<object.size(); i++) {
        Particle& p = object[i];
        Vector2f pic, flip = p.velocity;
        Matrix2f& grad = p.velocity_gradient;
        //grad.setData(0.0);
        grad = grad.setZero();
        float ox = p.grid_position(0), oy = p.grid_position(1);
        for (int idx=0, y=oy-1, y_end=y+3; y<=y_end; y++) {
            for (int x=ox-1, x_end=x+3; x<=x_end; x++, idx++) {
                float w = p.weights[idx];
                if (w > BSPLINE_EPSILON) {
                    GridNode &node = nodes[(int) (y*size(0)+x)];
                    pic += node.force*w;
                    flip += (node.force - node.velocity)*w;
                    grad += node.force * p.weight_gradient[idx].transpose();
                }
            }
        }
        p.velocity = flip*FLIP_PERCENT + pic*(1-FLIP_PERCENT);
    }
    collisionParticles();
}

void Grid::collisionGrid() {
    Vector2f delta_scale = Vector2f(timestep, timestep);
    delta_scale(0) /= cellsize(0);
    delta_scale(1) /= cellsize(1);
    for (int y=0, idx=0; y<size(1); y++) {
        for (int x=0; x<size(0); x++, idx++) {
            GridNode &node = nodes[idx];
            if (node.on) {
                Vector2f scaled_force;
                scaled_force(0) = node.force(0)*delta_scale(0);
                scaled_force(1) = node.force(1)*delta_scale(1);
                Vector2f new_pos = scaled_force+Vector2f(x, y);
                // LR
                if (new_pos(0) < BSPLINE_RADIUS || new_pos(0) > size(0)-BSPLINE_RADIUS-1) {
                    node.force(0) = 0;
                    node.force(1) *= STICKY;
                }
                // TB
                if (new_pos(1) < BSPLINE_RADIUS || new_pos(1) > size(1)-BSPLINE_RADIUS-1) {
                    node.force(1) = 0;
                    node.force(0) *= STICKY;
                }
            }
        }
    }
}

void Grid::collisionParticles() {
    for (int i = 0; i < object.size(); i++) {
        Particle &p = object[i];
        Vector2f temp_vel = Vector2f(p.velocity(0) / cellsize(0), p.velocity(1) / cellsize(1));
        Vector2f new_pos = p.grid_position + timestep * temp_vel;
        // LR
        if (new_pos(0) < BSPLINE_RADIUS - 1 || new_pos(0) > size(0) - BSPLINE_RADIUS)
            p.velocity(0) = -STICKY * p.velocity(0);
        // TB
        if (new_pos(1) < BSPLINE_RADIUS - 1 || new_pos(1) > size(1) - BSPLINE_RADIUS)
            p.velocity(1) = -STICKY * p.velocity(1);

    }
}