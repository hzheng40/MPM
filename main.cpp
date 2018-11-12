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
    Poisson poisson_sampler = Poisson(200, 200, 2.0, 30);
    poisson_sampler.initGrid();
    poisson_sampler.sample();
    poisson_sampler.writePartio("poisson.bgeo");
	return 0;
}