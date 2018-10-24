//
// Created by billyzheng on 10/23/18.
//

#ifndef MPM_H
#define MPM_H

#endif //MPM_H

#include <Partio.h>
#include <Eigen/Dense>
#include <iostream>

namespace MPM {
    class Particle;
    class ParticleList;
    class Grid;

    class Grid {
    public:
        Grid(int size, float resolution);
        ~Grid();
        void updateGrid(Particle& particle);

    private:
        int grid_size;
        float grid_resolution;
        Eigen::MatrixXi grid_ind;
        Eigen::MatrixXf grid_mass;
        Eigen::MatrixXf grid_vel;

    };

    class Particle {
    public:
        Particle();
        ~Particle();
        void initParticle();
        void updateParticle();
        void writeParticlesToFile(std::string file_path);

    private:
        float mass;
        float velocities;
        float volume;
        // 2d position, Vector3f if 3D
        Eigen::Vector2f position;
        Eigen::MatrixXf deformation_gradient;
    };

    class ParticleList {
    public:
        ParticleList();
        ~ParticleList();
        void init();

    private:

    };

//    float w interpolate()
}