//
// Created by billyzheng on 10/23/18.
//

#ifndef MPM_H
#define MPM_H

#endif //MPM_H

#include <Partio.h>

namespace MPM {
    class Particle;
    class Grid;

    class Grid {
    public:
        Grid(int size, float resolution);
        ~Grid();
        void initGrid();
        void updateGrid(Particles particles);

    private:
        int gird_size;
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
        void updateParticle(Grid grid);

    private:
        float mass;
        float velocities;
    };
}