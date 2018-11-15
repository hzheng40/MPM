//
// Created by billyzheng on 10/24/18.
//

#include <iostream>
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
    timestep = TIMESTEP;
}

Grid::~Grid() {
    delete[] nodes;
}

// map particle mass to grid
void Grid::initializeMass() {
    for (int i=0; i<nodes_length; i++) {
        GridNode &node = nodes[i];
        node.mass = 0.0;
    }
    for (int i=0; i<object.size(); i++) {
        Particle &p = object[i];
        Vector3f particle_dist = p.position - origin;
        p.grid_position(0) = particle_dist(0) / cellsize(0);
        p.grid_position(1) = particle_dist(1) / cellsize(1);
        p.grid_position(2) = particle_dist(2) / cellsize(2);
        int px_grid = static_cast<int>(p.grid_position(0)),
                py_grid = static_cast<int>(p.grid_position(1)),
                pz_grid = static_cast<int>(p.grid_position(2));
        for (int idx = 0, x_offset = -1; x_offset <= 2; x_offset++) {
            // x: column
            float g_node_x = (px_grid + x_offset) - p.grid_position(0);
            float wx = Grid::bspline(g_node_x), dx = Grid::bsplinePrime(g_node_x);
            for (int y_offset = -1; y_offset <= 2; y_offset++) {
                // y: row
                float g_node_y = (py_grid + y_offset) - p.grid_position(1);
                float wy = Grid::bspline(g_node_y), dy = Grid::bsplinePrime(g_node_y);
                for (int z_offset = -1; z_offset <= 2; z_offset++, idx++) {
                    // z: depth
                    float g_node_z = (pz_grid + z_offset) - p.grid_position(2);
                    float wz = Grid::bspline(g_node_z), dz = Grid::bsplinePrime(g_node_z);
                    // precomputing weight & weight gradient
                    float weight = wx * wy * wz;
                    p.weights[idx] = weight;
                    p.weight_gradient[idx] = Vector3f((dx * wy * wz) / cellsize(0),
                                                      (wx * dy * wz) / cellsize(1),
                                                      (wx * wy * dz) / cellsize(2));
                    // compute grid mass
                    int x = px_grid + x_offset, y = py_grid + y_offset, z = pz_grid + z_offset;
                    if (x<0 || x>size(0)-1 || y<0 || y>size(1)-1 || z<0 || z>size(2)-1) {
                        // TODO: prolly not correct to continue
                        continue;
                    }
                    nodes[(int) (y * size(0) + x + z * size(0) * size(1))].mass += weight * p.mass;
                }
            }
        }
    }
}

// map particle velocity to grid
void Grid::initializeVel() {
    for (int i=0; i<object.size(); i++) {
        Particle& p = object[i];
        Vector3f particle_dist = p.position - origin;
        p.grid_position(0) = particle_dist(0)/cellsize(0);
        p.grid_position(1) = particle_dist(1)/cellsize(1);
        p.grid_position(2) = particle_dist(2)/cellsize(2);
        int px_grid = static_cast<int>(p.grid_position(0)),
            py_grid = static_cast<int>(p.grid_position(1)),
            pz_grid = static_cast<int>(p.grid_position(2));
        for (int idx = 0, x_offset=-1; x_offset <= 2; x_offset++) {
            for (int y_offset=-1; y_offset <= 2; y_offset++) {
                for (int z_offset=-1; z_offset <= 2; z_offset++, idx++) {
                    float w = p.weights[idx];
                    if (w > BSPLINE_EPSILON) {
                        int x = px_grid+x_offset, y = py_grid+y_offset, z = pz_grid+z_offset;
                        int ind = (int) (y * size(0) + x + z * size(0) * size(1));
                        nodes[ind].velocity += p.velocity * w * p.mass;
                        nodes[ind].on = true;
                    }
                }
            }
        }
    }
    for (int i = 0; i < nodes_length; i++) {
        GridNode &node = nodes[i];
        if (node.on) {
            node.velocity /= node.mass;
        }
    }
    collisionGrid();
}

// initializing grid mass from original particles, only ran once at the beginning
void Grid::initialize() {
//    memset(nodes, 0, sizeof(GridNode)*nodes_length);
    for (int i=0; i<object.size(); i++) {
        Particle& p = object[i];
        Vector3f particle_dist = p.position - origin;
        p.grid_position(0) = particle_dist(0)/cellsize(0);
        p.grid_position(1) = particle_dist(1)/cellsize(1);
        p.grid_position(2) = particle_dist(2)/cellsize(2);
        int px_grid = static_cast<int>(p.grid_position(0)),
            py_grid = static_cast<int>(p.grid_position(1)),
            pz_grid = static_cast<int>(p.grid_position(2));
        for (int idx = 0, x_offset=-1; x_offset <= 2; x_offset++) {
            // x: column
            float g_node_x  = (px_grid+x_offset) - p.grid_position(0);
            float wx = Grid::bspline(g_node_x), dx = Grid::bsplinePrime(g_node_x);
            for (int y_offset=-1; y_offset <= 2; y_offset++) {
                // y: row
                float g_node_y = (py_grid+y_offset) - p.grid_position(1);
                float wy = Grid::bspline(g_node_y), dy = Grid::bsplinePrime(g_node_y);
                for (int z_offset=-1; z_offset <= 2; z_offset++, idx++) {
                // z: depth
                    float g_node_z = (pz_grid+z_offset) - p.grid_position(2);
                    float wz = Grid::bspline(g_node_z), dz = Grid::bsplinePrime(g_node_z);
                    // precomputing weight & weight gradient
                    float weight = wx*wy*wz;
                    p.weights[idx] = weight;
                    p.weight_gradient[idx] = Vector3f((dx*wy*wz)/cellsize(0), (wx*dy*wz)/cellsize(1), (wx*wy*dz)/cellsize(2));
                    // compute grid mass
                    int x = px_grid+x_offset, y = py_grid+y_offset, z = pz_grid+z_offset;
                    nodes[(int) (y*size(0)+x+z*size(0)*size(1))].mass += weight*p.mass;
                    if (weight > BSPLINE_EPSILON) {
                        int ind = (int) (y * size(0) + x + z * size(0) * size(1));
                        nodes[ind].velocity += p.velocity * weight * p.mass;
                        nodes[ind].on = true;
                    }
                }
            }
        }
    }
    for (int i = 0; i < nodes_length; i++) {
        GridNode &node = nodes[i];
        if (node.on) {
            node.velocity /= node.mass;
        }
    }
    collisionGrid();
}

// volumn from grid to particles, ran once at the start
void Grid::calculateVolumes() {
    for (int i=0; i<object.size(); i++) {
        Particle& p = object[i];
        int px_grid = static_cast<int>(p.grid_position(0)),
            py_grid = static_cast<int>(p.grid_position(1)),
            pz_grid = static_cast<int>(p.grid_position(2));
        p.density = 0.0;
        for (int idx = 0, x_offset=-1; x_offset <= 2; x_offset++) {
            for (int y_offset = -1; y_offset <= 2; y_offset++) {
                for (int z_offset = -1; z_offset <= 2; z_offset++, idx++) {
                    int x = px_grid+x_offset, y = py_grid+y_offset, z = pz_grid+z_offset;
                    float w = p.weights[idx];
                    if (w > BSPLINE_EPSILON) {
                        p.density += w*nodes[static_cast<int>(y*size(0)+x+z*size(0)*size(1))].mass;
                    }
                }
            }
        }
        p.density /= node_volume;
        p.volume = p.mass/p.density;
    }
}
// p2g
void Grid::p2g_vel(const Vector3f &gravity) {
    for (int i = 0; i < object.size(); i++) {
        Particle &p = object[i];
        Matrix3f energy = p.energyDerivative();
        int px_grid = static_cast<int>(p.grid_position(0)),
                py_grid = static_cast<int>(p.grid_position(1)),
                pz_grid = static_cast<int>(p.grid_position(2));
        for (int idx = 0, x_offset = -1; x_offset <= 2; x_offset++) {
            for (int y_offset = -1; y_offset <= 2; y_offset++) {
                for (int z_offset = -1; z_offset <= 2; z_offset++, idx++) {
                    int x = px_grid + x_offset, y = py_grid + y_offset, z = pz_grid + z_offset;
                    float w = p.weights[idx];
                    if (w > BSPLINE_EPSILON) {
                        int ind = (int) (y * size(0) + x + z * size(0) * size(1));
                        nodes[ind].new_velocity += energy * p.weight_gradient[idx];
                    }
                }
            }
        }
    }

    for (int i = 0; i < nodes_length; i++) {
        GridNode &node = nodes[i];
        if (node.on) {
            node.new_velocity = node.velocity + timestep * (gravity - node.new_velocity / node.mass);
        }
    }
    collisionGrid();
}

// g2p
void Grid::g2p_vel() {
    for (int i = 0; i < object.size(); i++) {
        Particle &p = object[i];
        Vector3f pic, flip = p.velocity;
        Matrix3f &grad = p.velocity_gradient;
        grad = grad.setZero();
        int px_grid = static_cast<int>(p.grid_position(0)),
            py_grid = static_cast<int>(p.grid_position(1)),
            pz_grid = static_cast<int>(p.grid_position(2));
        for (int idx = 0, x_offset = -1; x_offset <= 2; x_offset++) {
            for (int y_offset = -1; y_offset <= 2; y_offset++) {
                for (int z_offset = -1; z_offset <= 2; z_offset++, idx++) {
                    int x = px_grid + x_offset, y = py_grid + y_offset, z = pz_grid + z_offset;
                    float w = p.weights[idx];
                    if (w > BSPLINE_EPSILON) {
                        int ind = (int) (y * size(0) + x + z * size(0) * size(1));
                        GridNode &node = nodes[ind];
                        pic += node.new_velocity*w;
                        flip += (node.new_velocity - node.velocity)*w;
                        // grad is outer product
                        grad += node.new_velocity * p.weight_gradient[idx].transpose();
                    }
                }
            }
        }
        p.velocity = flip*FLIP_PERCENT + pic*(1-FLIP_PERCENT);
//        float vel_x = p.velocity(0);
//        float vel_y = p.velocity(1);
//        float vel_z = p.velocity(2);
//        cout << "part vel: (" << vel_x << ", " << vel_y << ", " << vel_z << ") \n";
    }
    collisionParticles();
}


void Grid::collisionGrid() {
    Vector3f delta_scale = Vector3f(timestep, timestep, timestep);
    delta_scale(0) /= cellsize(0);
    delta_scale(1) /= cellsize(1);
    delta_scale(2) /= cellsize(2);
    for (int y=0, idx=0; y<size(1); y++) {
        for (int x=0; x<size(0); x++) {
            for (int z=0; z<size(2); z++, idx++) {
                GridNode &node = nodes[idx];
                if (node.on) {
                    Vector3f scaled_velocity;
                    scaled_velocity(0) = node.new_velocity(0) * delta_scale(0);
                    scaled_velocity(1) = node.new_velocity(1) * delta_scale(1);
                    scaled_velocity(2) = node.new_velocity(2) * delta_scale(2);
                    Vector3f new_pos = scaled_velocity + Vector3f(x, y, z);
                    // Left Right
                    if (new_pos(0) < BSPLINE_RADIUS || new_pos(0) > size(0) - BSPLINE_RADIUS - 1) {
                        node.new_velocity(0) = 0;
                        node.new_velocity(1) *= STICKY;
                        node.new_velocity(2) *= STICKY;
                    }
                    // Top Bottom
                    if (new_pos(1) < BSPLINE_RADIUS || new_pos(1) > size(1) - BSPLINE_RADIUS - 1) {
                        node.new_velocity(1) = 0;
                        node.new_velocity(0) *= STICKY;
                        node.new_velocity(2) *= STICKY;
                    }
                    // Front Back
                    if (new_pos(2) < BSPLINE_RADIUS || new_pos(2) > size(2) - BSPLINE_RADIUS - 1) {
                        node.new_velocity(2) = 0;
                        node.new_velocity(0) *= STICKY;
                        node.new_velocity(1) *= STICKY;
                    }
                }
            }
        }
    }
}

void Grid::collisionParticles() {
    for (int i = 0; i < object.size(); i++) {
        Particle &p = object[i];
        Vector3f temp_vel = Vector3f(p.velocity(0)/cellsize(0), p.velocity(1)/cellsize(1), p.velocity(2)/cellsize(2));
        Vector3f new_pos = p.grid_position + timestep * temp_vel;
        // Left Right
        if (new_pos(0) < BSPLINE_RADIUS - 1) {
            p.position(0) = BSPLINE_RADIUS;
            p.velocity(0) = 0;
            p.velocity(1) *= STICKY;
            p.velocity(2) *= STICKY;
        } else if (new_pos(0) > size(0) - BSPLINE_RADIUS) {
            p.position(0) = size(0) - BSPLINE_RADIUS;
            p.velocity(0) = 0;
            p.velocity(1) *= STICKY;
            p.velocity(2) *= STICKY;
        }
        // Top Bottom
        if (new_pos(1) < BSPLINE_RADIUS - 1) {
            p.position(1) = BSPLINE_RADIUS;
            p.velocity(1) = 0;
            p.velocity(0) *= STICKY;
            p.velocity(2) *= STICKY;
        } else if (new_pos(1) > size(1) - BSPLINE_RADIUS) {
            p.position(1) = size(1) - BSPLINE_RADIUS;
            p.velocity(1) = 0;
            p.velocity(0) *= STICKY;
            p.velocity(2) *= STICKY;
        }
        // Front Back
        if (new_pos(2) < BSPLINE_RADIUS - 1) {
            p.position(2) = BSPLINE_RADIUS;
            p.velocity(2) = 0;
            p.velocity(1) *= STICKY;
            p.velocity(0) *= STICKY;
        } else if (new_pos(2) > size(2) - BSPLINE_RADIUS) {
            p.position(2) = size(1) - BSPLINE_RADIUS;
            p.velocity(2) = 0;
            p.velocity(1) *= STICKY;
            p.velocity(0) *= STICKY;
        }
    }
}