//
// Created by billyzheng on 11/11/18.
//
#include "Poisson.h"
Poisson::Poisson(int grid_width, int grid_height, int grid_depth, float radius, int k) {
    for (int i=0; i<grid_width*grid_height*grid_depth; i++) grid.push_back(Vector3f(-1.0, -1.0, -1.0));
    this->grid_width = grid_width;
    this->grid_height = grid_height;
    this->grid_depth = grid_depth;
    this->radius = radius;
    this->k = k;
    cell_size = radius/sqrt(3.0);
    vector<int> active;
}

Poisson::~Poisson() {
}

void Poisson::initGrid() {
    // random device
    random_device rd;
    mt19937 mt(rd());
    uniform_real_distribution<double> x_dist(0.0, grid_width);
    uniform_real_distribution<double> y_dist(0.0, grid_height);
    uniform_real_distribution<double> z_dist(0.0, grid_depth);
    // random sample, x in [0, width), y in [0, height)
    float x = x_dist(mt);
    float y = y_dist(mt);
    float z = z_dist(mt);
    Vector3f init_pos = Vector3f(x, y, z);
    cout << "initial sample x: " << x << "\n";
    cout << "initial sample y: " << y << "\n";
    cout << "initial sample z: " << z << "\n";
    // figure out grid position of random point
    int x_i = static_cast<int> (y/cell_size);
    int y_i = static_cast<int> (x/cell_size);
    int z_i = static_cast<int> (z/cell_size);
    cout << "initial sample x grid pos: " << x_i << "\n";
    cout << "initial sample y grid pos: " << y_i << "\n";
    cout << "initial sample z grid pos: " << z_i << "\n";
    grid[x_i*grid_width + y_i + z_i*grid_width*grid_height] = init_pos;
    active.push_back(x_i*grid_width + y_i + z_i*grid_width*grid_height);
    cout << "done initialization" << "\n";
}

void Poisson::sample() {
    cout << "Starting sample" << "\n";
    // random device
    random_device rd;
    mt19937 mt(rd());

    uniform_real_distribution<double> theta_dist(0.0, 2*M_PI);
    uniform_real_distribution<double> phi_dist(0.0, M_PI);
    uniform_real_distribution<double> radius_dist(radius, 2*radius);

    // frame num for partio
    int frame_num = 0;

    while (!active.empty()) {
        uniform_int_distribution<int> idx_dist(0, active.size()-1);
//        cout << "loop num: " << frame_num << "\n";
//        cout << "active length: " << active.size() << "\n";
//        cout << "grid length: " << grid.size() << "\n";
//        cout << "active is empty: " << active.empty() << "\n";
        int rand_ind = idx_dist(mt);
        int active_ind = active[rand_ind];
        bool point_added = false;
        // k samples around current active point
        for (int i=0; i<k; i++) {
            // random point in donut around current active point
            float rand_theta = theta_dist(mt);
            float rand_phi = phi_dist(mt);
            float rand_r = radius_dist(mt);
            float rand_offset_x = rand_r*cos(rand_theta)*sin(rand_phi);
            float rand_offset_y = rand_r*sin(rand_theta)*sin(rand_phi);
            float rand_offset_z = rand_r*cos(rand_phi);
            Vector3f active_pt = grid[active_ind];
            float rand_pt_x = active_pt(0) + rand_offset_x;
            float rand_pt_y = active_pt(1) + rand_offset_y;
            float rand_pt_z = active_pt(2) + rand_offset_z;
//            Vector2f rand_pt = active_pt + Vector2f(rand_offset_x, rand_offset_y);
            // random point's grid position
            int rand_pt_i = static_cast<int> (rand_pt_y/cell_size);
            int rand_pt_j = static_cast<int> (rand_pt_x/cell_size);
            int rand_pt_k = static_cast<int> (rand_pt_z/cell_size);
//            cout << "current point we're checking (" << rand_pt_i << ", " << rand_pt_j << ")\n";
            bool valid = true;
            if (!(rand_pt_i >= 1 && rand_pt_i < grid_height-1
                    && rand_pt_j >= 1 && rand_pt_j < grid_width-1
                    && rand_pt_k >= 1 && rand_pt_k < grid_depth-1)) {
                continue;
            }
            // checking only neighbors
            Vector3f curr_pt_in_grid = grid[rand_pt_i*grid_width + rand_pt_j + rand_pt_k*grid_width*grid_height];
            if (curr_pt_in_grid(0) != -1.0 && curr_pt_in_grid(1) != -1.0 && curr_pt_in_grid(2) != -1.0) {
                continue;
            }
            for (int i_pt = -1; i_pt <= 1; i_pt++) {
                for (int j_pt = -1; j_pt <= 1; j_pt++) {
                    for (int k_pt = -1; k_pt <= 1; k_pt++) {
                        Vector3f neighbor = grid[(rand_pt_i+i_pt)*grid_width+(rand_pt_j+j_pt)+(rand_pt_k+k_pt)*grid_width*grid_height];
    //                    cout << "neighbor pos: (" << neighbor(0) << ", " << neighbor(1) << ")\n";
                        if (neighbor(0) != -1.0 && neighbor(1) != -1.0 && neighbor(2) != -1.0) {
                            float distance = sqrt(pow((neighbor(0) - rand_pt_x), 2)
                                    + pow((neighbor(1) - rand_pt_y), 2)
                                    + pow((neighbor(2) - rand_pt_z), 2));
                            if (distance < radius) valid = false;
                        }
                    }
                }
            }
            // if no violation add to grid and active list
            if (valid) {
                Vector3f point_to_add = Vector3f(rand_pt_x, rand_pt_y, rand_pt_z);
                point_added = true;
                grid[rand_pt_i * grid_width + rand_pt_j + rand_pt_k * grid_width * grid_height] = point_to_add;
                int new_active_ind = rand_pt_i * grid_width + rand_pt_j + rand_pt_k * grid_width * grid_height;
                active.push_back(new_active_ind);
            }
        }
        // remove from active list if no new point added
        if (!point_added) {
            if (!active.empty()) {
                active.erase(active.begin() + rand_ind );
            }
        }
        frame_num++;
//        writePartioByFrame("poisson", frame_num);
    }
    cout << "total loops ran: " << frame_num << "\n";
}

void Poisson::writePartio(const string& particle_file) {
    Partio::ParticlesDataMutable *parts = Partio::create();
    Partio::ParticleAttribute posH;
    posH = parts->addAttribute("position", Partio::VECTOR, 3);
    for (int i=0; i<grid_width*grid_height*grid_depth; i++) {
        Vector3f pt = grid[i];
        if (pt(0) >= 0 && pt(1) >= 0 && pt(2) >= 0) {
            int idx = parts->addParticle();
            float *p = parts->dataWrite<float>(posH, idx);
            p[0] = pt(0);
            p[1] = pt(1);
            p[2] = pt(2);
        }
    }
    Partio::write(particle_file.c_str(), *parts);
    parts->release();
}

void Poisson::writePartioByObejct(const string &particle_file,
                                  vector<Vector3f, aligned_allocator<Vector3f>> object_file){
    Partio::ParticlesDataMutable *parts = Partio::create();
    Partio::ParticleAttribute posH;
    posH = parts->addAttribute("position", Partio::VECTOR, 3);
    for (int i=0; i<object_file.size(); i++) {
        Vector3f pt = object_file[i];
        if (pt(0) >= 0 && pt(1) >= 0 && pt(2) >= 0) {
            int idx = parts->addParticle();
            float *p = parts->dataWrite<float>(posH, idx);
            p[0] = pt(0);
            p[1] = pt(1);
            p[2] = pt(2);
        }
    }
    Partio::write(particle_file.c_str(), *parts);
    parts->release();
}

void Poisson::writePartioByFrame(const string& particle_file_prefix, int seq_num) {
    Partio::ParticlesDataMutable *parts = Partio::create();
    Partio::ParticleAttribute posH;
    posH = parts->addAttribute("position", Partio::VECTOR, 3);
    for (int i=0; i<grid_width*grid_height*grid_depth; i++) {
        Vector3f pt = grid[i];
        if (pt(0) >= 0 && pt(1) >= 0 && pt(2) >= 0) {
            int idx = parts->addParticle();
            float *p = parts->dataWrite<float>(posH, idx);
            p[0] = pt(0);
            p[1] = pt(1);
            p[2] = pt(2);
        }
    }
    Partio::write((particle_file_prefix + to_string(seq_num) + ".bgeo").c_str(), *parts);
    parts->release();
}

vector<Vector3f, aligned_allocator<Vector3f>> Poisson::toSphere() {
    vector<Vector3f, aligned_allocator<Vector3f>> sphere_grid;
    float sphere_r = min(min(grid_width, grid_height), grid_depth) / 2.0;
    Vector3f sphere_c = Vector3f(grid_width/2.0, grid_height/2.0, grid_depth/2.0);
    for (int i=0; i<grid_width*grid_height*grid_depth; i++) {
        Vector3f curr_pt = grid[i];
        float dist = pow((curr_pt(0)-sphere_c(0)), 2) + pow((curr_pt(1)-sphere_c(1)), 2) + pow((curr_pt(2)-sphere_c(2)), 2);
        if (dist <= pow(sphere_r, 2)) {
            sphere_grid.push_back(curr_pt);
        }
    }
    return sphere_grid;
}
vector<Particle> Poisson::toObject(vector<Vector3f, aligned_allocator<Vector3f>> object) {
    vector<Particle> part_list;
    for (int i=0; i<object.size(); i++) {
        Vector3f point = object[i];
        // create a particle for each point in pointcloud
        if (point(0) >= 0 && point(0) < grid_width
            && point(1) >= 0 && point(1) < grid_height
            && point(2) >= 0 && point(2) < grid_depth) {
            Particle particle = Particle(point, Vector3f::Zero(), PT_MASS, LAMBDA, MU, TIMESTEP);
            part_list.push_back(particle);
        } else {
            cout << "skipping \n";
        }
    }
    return part_list;
}