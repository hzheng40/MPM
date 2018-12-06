#include <Partio.h>
#include <Eigen/Dense>
#include <iostream>
#include "Constants.h"
#include "Particle.h"
#include "Grid.h"
#include "Poisson.h"

void writePartio(const string& particle_file_prefix, int seq_num, vector<Particle> part_list) {
    Partio::ParticlesDataMutable *parts = Partio::create();
    Partio::ParticleAttribute posH;
    posH = parts->addAttribute("position", Partio::VECTOR, 3);
    for (int i=0; i<part_list.size(); i++) {
        Particle &pt = part_list[i];
        int idx = parts->addParticle();
        float *p = parts->dataWrite<float>(posH, idx);
        p[0] = pt.position(0);
        p[1] = pt.position(1);
        p[2] = pt.position(2);
    }
    Partio::write((particle_file_prefix + to_string(seq_num) + ".bgeo").c_str(), *parts);
    parts->release();
}

int main() {

    Poisson poisson_sampler = Poisson(1, 1, 1, 0.03, 30);
    poisson_sampler.initGrid();
    poisson_sampler.sample();
    cout << "done sampling" << "\n";
//    poisson_sampler.writePartio("poisson.bgeo");
//    auto sphere = poisson_sampler.toSphere();
    auto cube = poisson_sampler.toCube();
//    poisson_sampler.writePartioByObejct("poisson_sphere.bgeo", sphere);
    poisson_sampler.writePartioByObejct("poisson_cube.bgeo", cube);
//    vector<Particle> part_list = poisson_sampler.toObject(sphere);
    vector<Particle> part_list = poisson_sampler.toObject(cube);
    cout << "particle in list now" << "\n";
    // make grid
    Grid grid(Vector3f(-0.5, -0.5, -0.5), Vector3f(2.0, 2.0, 2.0), Vector3f(100.0,100.0,100.0), part_list);
    grid.initializeMass();
    grid.calculateVolumes();
    Vector3f gravity = Vector3f(0, 0, GRAVITY);
    // main MPM loop
    for (int frame_num=0; frame_num<MAX_ITER; frame_num++) {
        cout << "current frame number: " << frame_num << "\n";
        grid.initializeMass();
        grid.initializeVel();
        grid.p2g_vel(gravity);
        for (int i=0; i<grid.object.size(); i++) {
            Particle &part = grid.object[i];
            part.updateGradient();
//            part.checksum();
            part.applyPlasticity();
        }
        grid.g2p_vel();
        for (int i=0; i<grid.object.size(); i++) {
            Particle &part = grid.object[i];
            part.updatePos();
        }
        writePartio("mpm_partio_3/mpm_", frame_num, grid.object);
    }
	return 0;
}