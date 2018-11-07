////
//// Created by Hongrui Zheng on 11/5/18.
////
#include <Partio.h>
#include <vector>
#include <string>
#include <cstdlib>
#include <ctime>
using T = float;
constexpr int dim = 3;

template <class T, int dim>
void writePartio(const std::string& particleFile)
{
    srand(static_cast<unsigned> (time(0)));
    std::array<float, 100> pointsx{};
    std::array<float, 100> pointsy{};
    std::array<float, 100> pointsz{};
    for (int idx=0; idx< 100; idx++) {
        pointsx[idx] = static_cast<float> (rand()) / (static_cast<float> (RAND_MAX/100.0));
        pointsy[idx] = static_cast<float> (rand()) / (static_cast<float> (RAND_MAX/100.0));
        pointsz[idx] = static_cast<float> (rand()) / (static_cast<float> (RAND_MAX/100.0));
    }
    for (int time=0; time<2000; time++) {
        std::array<float, 100> rands{};
        for (int j=0; j<100; j++) {
            pointsx[j] += static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
            pointsy[j] += static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
            pointsz[j] += static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

        }
        Partio::ParticlesDataMutable *parts = Partio::create();
        Partio::ParticleAttribute posH, vH, mH;
        mH = parts->addAttribute("m", Partio::VECTOR, 1);
        posH = parts->addAttribute("position", Partio::VECTOR, 3);
        vH = parts->addAttribute("v", Partio::VECTOR, 3);

        for (int i = 0; i < 100; i++) {
            int idx = parts->addParticle();
            float *m = parts->dataWrite<float>(mH, idx);
            float *p = parts->dataWrite<float>(posH, idx);
            float *v = parts->dataWrite<float>(vH, idx);
            m[0] = (T) i;
            p[0] = pointsx[i];
            p[1] = pointsy[i];
            p[2] = pointsz[i];
            for (int k = 0; k < 3; k++)
                v[k] = (T) i;
        }

        Partio::write((particleFile+std::to_string(time)+".bgeo").c_str(), *parts);
        parts->release();
    }
}

int main() {
    std::string file="test";
    writePartio<T,dim>(file);
    return 0;
}