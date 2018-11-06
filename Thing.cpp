//
// Created by billyzheng on 11/4/18.
//

#include "Thing.h"

//constructors


void Thing::toFile(std::string filename) {
    Partio::ParticlesDataMutable* parts = Partio::create();
    Partio::ParticleAttribute posH, vH, mH, fH, cH;
    mH = parts->addAttribute("m", Partio::VECTOR, 1);
    posH = parts->addAttribute("position", Partio::VECTOR, 3);
    vH = parts->addAttribute("v", Partio::VECTOR, 3);
    fH = parts->addAttribute("f", Partio::VECTOR, 3);
    cH = parts->addAttribute("color", Partio::VECTOR, 3);
    for (unsigned i=0; i< xs.size(); i++){
        int idx = parts->addParticle();
        float* m = parts->dataWrite<float>(mH, idx);
        float* p = parts->dataWrite<float>(posH, idx);
        float* v = parts->dataWrite<float>(vH, idx);
        float* f = parts->dataWrite<float>(fH, idx);
        float* c = parts->dataWrite<float>(cH, idx);
        m[0] = ms[i];
        for (int k = 0; k < 3; k++)
            p[k] = (T)xs[i][k];
        for (int k = 0; k < 3; k++)
            v[k] = (T)vs[i][k];
        for (int k = 0; k < 3; k++)
            f[k] = (T)fs[i][k] / 10.f;
        //c[k] = (T)colors[i][k];
    }

    Partio::write(filename.c_str(), *parts);
    parts->release();
}