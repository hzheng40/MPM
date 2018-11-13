#include <Partio.h>
#include <Eigen/Dense>
#include "Constants.h"
//#include "Particle.h"
//#include "Grid.h"
#include "Poisson.h"

//Grid initGrid() {
//    Grid grid = Grid::Grid()
//}

int main() {
//    Grid grid = initGrid();
//    for (int time_step=0; time_step<MAX_TIMESTEP; time_step++){
//
//    }
    Poisson poisson_sampler = Poisson(100, 100, 100, 2.0, 30);
    poisson_sampler.initGrid();
    poisson_sampler.sample();
    cout << "done sampling" << "\n";
    poisson_sampler.writePartio("poisson.bgeo");
    auto sphere = poisson_sampler.toSphere();
    poisson_sampler.writePartioByObejct("poisson_sphere.bgeo", sphere);
	return 0;
}