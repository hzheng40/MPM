//
// Created by billyzheng on 10/24/18.
//

#ifndef MPM_GRID_H
#define MPM_GRID_H

#include <Eigen/Dense>
#include "Constants.h"

namespace MPM {
    typedef struct GridNode {
        float mass;
        Eigen::Vector2f velocity;

    } GridNode;

    class Grid {
    public:
    };

} // end namespace

#endif //MPM_GRID_H
