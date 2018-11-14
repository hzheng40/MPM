#include <Partio.h>
#include <Eigen/Dense>
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

    Poisson poisson_sampler = Poisson(50, 50, 50, 2.0, 30);
    poisson_sampler.initGrid();
    poisson_sampler.sample();
    cout << "done sampling" << "\n";
//    poisson_sampler.writePartio("poisson.bgeo");
    auto sphere = poisson_sampler.toSphere();
//    poisson_sampler.writePartioByObejct("poisson_sphere.bgeo", sphere);
    vector<Particle> part_list = poisson_sampler.toObject(sphere);
    cout << "particle in list now" << "\n";
    // make grid
    Grid grid(Vector3f(-100.0,-100.0,-100.0), Vector3f(200.0,200.0,200.0), Vector3f(200.0,200.0,200.0), part_list);
    grid.initialize();
    grid.calculateVolumes();
    Vector3f gravity = Vector3f(0, 0, GRAVITY);
    // main MPM loop
    for (float time_step=0; time_step<MAX_TIMESTEP; time_step+=TIMESTEP) {
        int frame_num = static_cast<int> (time_step / TIMESTEP);
        grid.initialize();
        grid.p2g_vel(gravity);
        grid.g2p_vel();
        for (int i=0; i<part_list.size(); i++) {
            Particle &part = part_list[i];
            part.updatePos();
            part.updateGradient();
            part.applyPlasticity();
        }
        writePartio("mpm_", frame_num, part_list);
    }
	return 0;
}