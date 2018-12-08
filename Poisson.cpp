//
// Created by billyzheng on 11/11/18.
//
#include "Poisson.h"
Poisson::Poisson(double grid_width, double grid_height, double grid_depth, double radius, int k) {
    cell_size = radius/sqrt(3.0);
    this->grid_width = grid_width;
    this->grid_height = grid_height;
    this->grid_depth = grid_depth;
    this->radius = radius;
    this->k = k;
    grid_length = (int) (grid_width/cell_size)*(grid_height/cell_size)*(grid_depth/cell_size);
    w_ind = grid_width/cell_size;
    h_ind = grid_height/cell_size;
    d_ind = grid_depth/cell_size;
    for (int i=0; i<grid_length; i++) grid.push_back(Vector3d(-1.0, -1.0, -1.0));
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
    double x = x_dist(mt);
    double y = y_dist(mt);
    double z = z_dist(mt);
    Vector3d init_pos = Vector3d(x, y, z);
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

    grid[x_i*w_ind + y_i + z_i*w_ind*h_ind] = init_pos;
    active.push_back(x_i*w_ind + y_i + z_i*w_ind*h_ind);
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
            double rand_theta = theta_dist(mt);
            double rand_phi = phi_dist(mt);
            double rand_r = radius_dist(mt);
            double rand_offset_x = rand_r*cos(rand_theta)*sin(rand_phi);
            double rand_offset_y = rand_r*sin(rand_theta)*sin(rand_phi);
            double rand_offset_z = rand_r*cos(rand_phi);
            Vector3d active_pt = grid[active_ind];
            double rand_pt_x = active_pt(0) + rand_offset_x;
            double rand_pt_y = active_pt(1) + rand_offset_y;
            double rand_pt_z = active_pt(2) + rand_offset_z;
//            Vector2f rand_pt = active_pt + Vector2f(rand_offset_x, rand_offset_y);
            // random point's grid position
            int rand_pt_i = static_cast<int> (rand_pt_y/cell_size);
            int rand_pt_j = static_cast<int> (rand_pt_x/cell_size);
            int rand_pt_k = static_cast<int> (rand_pt_z/cell_size);
//            cout << "current point we're checking (" << rand_pt_i << ", " << rand_pt_j << ")\n";
            bool valid = true;
            if (!(rand_pt_i >= 1 && rand_pt_i < h_ind-1
                    && rand_pt_j >= 1 && rand_pt_j < w_ind-1
                    && rand_pt_k >= 1 && rand_pt_k < d_ind-1)) {
                continue;
            }
            // checking only neighbors
            Vector3d curr_pt_in_grid = grid[rand_pt_i*w_ind + rand_pt_j + rand_pt_k*w_ind*h_ind];
            if (curr_pt_in_grid(0) != -1.0 && curr_pt_in_grid(1) != -1.0 && curr_pt_in_grid(2) != -1.0) {
                continue;
            }
            for (int i_pt = -1; i_pt <= 1; i_pt++) {
                for (int j_pt = -1; j_pt <= 1; j_pt++) {
                    for (int k_pt = -1; k_pt <= 1; k_pt++) {
                        Vector3d neighbor = grid[(rand_pt_i+i_pt)*w_ind+(rand_pt_j+j_pt)+(rand_pt_k+k_pt)*w_ind*h_ind];
    //                    cout << "neighbor pos: (" << neighbor(0) << ", " << neighbor(1) << ")\n";
                        if (neighbor(0) != -1.0 && neighbor(1) != -1.0 && neighbor(2) != -1.0) {
                            double distance = sqrt(pow((neighbor(0) - rand_pt_x), 2)
                                    + pow((neighbor(1) - rand_pt_y), 2)
                                    + pow((neighbor(2) - rand_pt_z), 2));
                            if (distance < radius) valid = false;
                        }
                    }
                }
            }
            // if no violation add to grid and active list
            if (valid) {
                Vector3d point_to_add = Vector3d(rand_pt_x, rand_pt_y, rand_pt_z);
                point_added = true;
                grid[rand_pt_i * w_ind + rand_pt_j + rand_pt_k * w_ind * h_ind] = point_to_add;
                int new_active_ind = rand_pt_i * w_ind + rand_pt_j + rand_pt_k * w_ind * h_ind;
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
        Vector3d pt = grid[i];
        if (pt(0) >= 0 && pt(1) >= 0 && pt(2) >= 0) {
            int idx = parts->addParticle();
            double *p = parts->dataWrite<double>(posH, idx);
            p[0] = pt(0);
            p[1] = pt(1);
            p[2] = pt(2);
        }
    }
    Partio::write(particle_file.c_str(), *parts);
    parts->release();
}

void Poisson::writePartioByObejct(const string &particle_file,
                                  vector<Vector3d, aligned_allocator<Vector3d>> object_file){
    Partio::ParticlesDataMutable *parts = Partio::create();
    Partio::ParticleAttribute posH;
    posH = parts->addAttribute("position", Partio::VECTOR, 3);
    for (int i=0; i<object_file.size(); i++) {
        Vector3d pt = object_file[i];
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
        Vector3d pt = grid[i];
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

vector<Vector3d, aligned_allocator<Vector3d>> Poisson::toSphere() {
    vector<Vector3d, aligned_allocator<Vector3d>> sphere_grid;
    double sphere_r = min(min(grid_width, grid_height), grid_depth) / 2.0;
    Vector3d sphere_c = Vector3d(grid_width/2.0, grid_height/2.0, grid_depth/2.0);
    for (int i=0; i<w_ind*h_ind*d_ind; i++) {
        Vector3d curr_pt = grid[i];
        double dist = pow((curr_pt(0)-sphere_c(0)), 2) + pow((curr_pt(1)-sphere_c(1)), 2) + pow((curr_pt(2)-sphere_c(2)), 2);
        if (dist <= pow(sphere_r, 2)) {
            sphere_grid.push_back(curr_pt);
        }
    }
    return sphere_grid;
}

vector<Vector3d, aligned_allocator<Vector3d>> Poisson::toCube() {
    vector<Vector3d, aligned_allocator<Vector3d>> sphere_grid;
    for (int i=0; i<w_ind*h_ind*d_ind; i++) {
        Vector3d curr_pt = grid[i];
        sphere_grid.push_back(curr_pt);
    }
    return sphere_grid;
}

vector<Particle> Poisson::toObject(vector<Vector3d, aligned_allocator<Vector3d>> object) {
    vector<Particle> part_list;
    for (int i=0; i<object.size(); i++) {
        Vector3d point = object[i];
        // create a particle for each point in pointcloud
        if (point(0) >= 0 && point(0) < grid_width
            && point(1) >= 0 && point(1) < grid_height
            && point(2) >= 0 && point(2) < grid_depth) {
//            Particle particle = Particle(point, Vector3d::Zero(), PT_MASS, LAMBDA, MU, TIMESTEP);
            Particle particle = Particle(point, Vector3d(0,0,0), PT_MASS, LAMBDA, MU, TIMESTEP);
            part_list.push_back(particle);
        } else {
//            cout << "skipping \n";
            continue;
        }
    }
    return part_list;
}