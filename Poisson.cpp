//
// Created by billyzheng on 11/11/18.
//

#include "Poisson.h"
#include "Partio.h"
#include <random>
#include <iostream>

Poisson::Poisson(int grid_width, int grid_height, float radius, int k) {
    for (int i=0; i<grid_width*grid_height; i++) grid.emplace_back(Vector2f(-1, -1));
    this->grid_width = grid_width;
    this->grid_height = grid_height;
    this->radius = radius;
    this->k = k;
    cell_size = radius/sqrt(2.0);

}

Poisson::~Poisson() {
}

void Poisson::initGrid() {
    // random device
    random_device rd;
    mt19937 mt(rd());
    uniform_real_distribution<double> x_dist(0.0, grid_width);
    uniform_real_distribution<double> y_dist(0.0, grid_height);
    // random sample, x in [0, width), y in [0, height)
    float x = x_dist(mt);
    float y = y_dist(mt);
    Vector2f init_pos = Vector2f(x, y);
    cout << "initial sample x: " << x << "\n";
    cout << "initial sample y: " << y << "\n";
    // figure out grid position of random point
    int x_i = static_cast<int> (x/cell_size);
    int y_i = static_cast<int> (y/cell_size);
    cout << "initial sample x grid pos: " << x_i << "\n";
    cout << "initial sample y grid pos: " << y_i << "\n";
    grid[x_i+y_i*grid_width] = init_pos;
    active.push_back(init_pos);
    cout << "done initialization" << "\n";
}

void Poisson::sample() {
    cout << "Starting sample" << "\n";
    // random device
    random_device rd;
    mt19937 mt(rd());
    uniform_int_distribution<int> idx_dist(0, active.size());
    uniform_real_distribution<double> ang_dist(0.0, 2*M_PI);
    uniform_real_distribution<double> radius_dist(radius, 2*radius);

    // frame num for partio
    int frame_num = 0;

    while (!active.empty()) {
        cout << "loop num: " << frame_num << "\n";
        int rand_ind = idx_dist(mt);
        Vector2f active_pt = active[rand_ind];
        bool point_added = false;
        // k samples around current active point
        for (int i=0; i<k; i++) {
            // random point in donut around current active point
            float rand_angle = ang_dist(mt);
            float rand_r = radius_dist(mt);
            float rand_pt_x = rand_r*cos(rand_angle);
            float rand_pt_y = rand_r*sin(rand_angle);
            Vector2f rand_pt = active_pt + Vector2f(rand_pt_x, rand_pt_y);
            // random point's grid position
            int rand_pt_i = static_cast<int> (rand_pt(0)/cell_size);
            int rand_pt_j = static_cast<int> (rand_pt(1)/cell_size);
            bool valid = true;
            // checking only neighbors
            for (int i_pt = -1; i_pt <= 1; i_pt++) {
                for (int j_pt = -1; j_pt <= 1; j_pt++) {
                    Vector2f neighbor = grid[(rand_pt_i+i_pt)+(rand_pt_j+j_pt)*grid_width];
                    if (neighbor(0) != -1.0 && neighbor(1) != -1.0) {
                        float distance = sqrt(pow((neighbor(0) - rand_pt(0)), 2) + pow((neighbor(1) - rand_pt(1)), 2));
                        if (distance < radius) valid = false;
                    }
                }
            }
            // if no violation add to grid and active list
            if (valid) {
                point_added = true;
                active.push_back(rand_pt);
                grid[rand_pt_i + rand_pt_j * grid_width] = rand_pt;
            }
        }
        // remove from active list if no new point added
        if (!point_added) {
            active.erase(active.begin()+rand_ind-1);
        }
        frame_num++;

    }
}

void Poisson::writePartio(const string& particle_file) {
    Partio::ParticlesDataMutable *parts = Partio::create();
    Partio::ParticleAttribute posH;
    posH = parts->addAttribute("position", Partio::VECTOR, 3);
    for (int i=0; i<grid_width*grid_height; i++) {
        int idx = parts->addParticle();
        float *p = parts->dataWrite<float>(posH, idx);
        Vector2f pt = grid[i];
        if (pt(0) != -1.0 && pt(1) != -1.0) {
            p[0] = pt(0);
            p[1] = pt(1);
            p[2] = 0.0;
        }
    }
    Partio::write(particle_file.c_str(), *parts);
    parts->release();
}

//void Poisson::toObject(array& object) {
//
//}