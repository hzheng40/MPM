////
//// Created by Hongrui Zheng on 11/5/18.
////
//#include <Partio.h>
//#include <Eigen/Dense>
//using namespace Eigen;
//using namespace std;
//
//int main() {
//    string filename = "rand";
//    Partio::ParticlesDataMutable* parts = Partio::create();
//    Partio::ParticleAttribute posH, vH, mH, fH, cH;
//    mH = parts->addAttribute("m", Partio::VECTOR, 1);
//    posH = parts->addAttribute("position", Partio::VECTOR, 3);
//    vH = parts->addAttribute("v", Partio::VECTOR, 3);
//    fH = parts->addAttribute("f", Partio::VECTOR, 3);
//    cH = parts->addAttribute("color", Partio::VECTOR, 3);
//    for (unsigned i=0; i<10; i++) {
//        for (unsigned j = 0; j < 10; j++) {
//            int idx = parts->addParticle();
//            float *m = parts->dataWrite<float>(mH, idx);
//            float *p = parts->dataWrite<float>(posH, idx);
//            float *v = parts->dataWrite<float>(vH, idx);
//            float *f = parts->dataWrite<float>(fH, idx);
//            float *c = parts->dataWrite<float>(cH, idx);
//            m[0] = 100.0;
//            p[0] = i; p[1] = j;
//
////            for (int k = 0; k < 3; k++)
////                p[k] = (T) xs[i][k];
////            for (int k = 0; k < 3; k++)
////                v[k] = (T) vs[i][k];
////            for (int k = 0; k < 3; k++)
////                f[k] = (T) fs[i][k] / 10.f;
//            //c[k] = (T)colors[i][k];
//        }
//    }
//    Partio::write(filename.c_str(), *parts);
//    parts->release();
//}


#include <Partio.h>
#include <vector>
#include <string>

using T = float;
constexpr int dim = 3;

template <class T, int dim>
void writePartio(const std::string& particleFile)
{
    Partio::ParticlesDataMutable* parts = Partio::create();
    Partio::ParticleAttribute posH, vH, mH;
    mH = parts->addAttribute("m", Partio::VECTOR, 1);
    posH = parts->addAttribute("position", Partio::VECTOR, 3);
    vH = parts->addAttribute("v", Partio::VECTOR, 3);
    for (int i=0; i<3; i++){
        int idx = parts->addParticle();
        float* m = parts->dataWrite<float>(mH, idx);
        float* p = parts->dataWrite<float>(posH, idx);
        float* v = parts->dataWrite<float>(vH, idx);
        m[0] = (T)i;
        for (int k = 0; k < 3; k++)
            p[k] = (T)i;
        for (int k = 0; k < 3; k++)
            v[k] = (T)i;
    }

    Partio::write(particleFile.c_str(), *parts);
    parts->release();
}

int main() {
    std::string file="test.bgeo";
    writePartio<T,dim>(file);
    return 0;
}