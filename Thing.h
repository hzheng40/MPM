//
// Created by billyzheng on 11/4/18.
//

#ifndef MPM_OBJECT_H
#define MPM_OBJECT_H
#include <vector>
#include <Eigen/Dense>
#include "Particle.h"
#include "Constants.h"
using namespace Eigen;
class Thing {
public:
    int size;
    float max_velocity;
    std::vector<Particle> particles;
    Thing();
    Thing(int object_size);
    Thing(const Thing& orig);
    virtual ~Thing();
    // transformations
    void scale(Vector2f origin, Vector2f scale);
    void translate(Vector2f offset);
    // update particle area
    void update();
    // get bounding box
    void bounds(float bound[4]);
    
};
#endif //MPM_OBJECT_H
