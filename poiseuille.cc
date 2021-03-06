#include <iostream>
#include <vector>
#include <algorithm>

#include "lbm.h"
#include "boundary_conditions.h"

constexpr std::size_t dimY = 20;
constexpr std::size_t dimX = 5*dimY;

constexpr double uInflow  = 0.02;
constexpr double reynolds = 100;

constexpr double tau   = 3. * uInflow * (dimX-1) / reynolds + 0.5;
constexpr double omega = 1. / tau;

DataCellBuffer pop(dimX, dimY);
FluidBuffer fluid(dimX, dimY);

double poiseuille(std::size_t y) {
	return -4. * uInflow / (dimY*dimY) * (y+0.5) * (y+0.5 - dimY);
}

void init() {
	for ( std::size_t x = 0; x < dimX; ++x ) {
		for ( std::size_t y = 0; y < dimY; ++y ) {
			fluid.density(x,y)  = 1.0;
			fluid.velocity(x,y) = {0.0, 0.0};

			static_cast<Cell&>(pop.curr(x,y)).equilibrize(
				fluid.density(x,y), fluid.velocity(x,y));
			static_cast<Cell&>(pop.prev(x,y)).equilibrize(
				fluid.density(x,y), fluid.velocity(x,y));
		}
	}
}

void computeLbmStep() {
	pop.swap();

	for ( std::size_t x = 0; x < dimX; ++x ) {
		for ( std::size_t y = 0; y < dimY; ++y ) {
			streamFluidCell(pop, x, y);
		}
	}

	for ( std::size_t x = 0; x < dimX; ++x ) {
		computeWallCell(pop, {x, 0     }, { 0, 1});
		computeWallCell(pop, {x, dimY-1}, { 0,-1});
	}

	for ( std::size_t y = 1; y < dimY-1; ++y ) {
		computeMovingWallCell(pop, {0,y}, {1,0}, {poiseuille(y), 0});
	}

	for ( std::size_t x = 0; x < dimX; ++x ) {
		for ( std::size_t y = 0; y < dimY; ++y ) {
			Cell& cell = static_cast<Cell&>(pop.curr(x,y));

			// bulk density
			fluid.density(x,y) = cell.sum();

			// outflow density condition
			if ( x == dimX-1 && y > 0 && y < dimY-1 ) {
				fluid.density(x,y) = 1.0;
			}

			fluid.velocity(x,y) = cell.velocity(fluid.density(x,y));

			collideFluidCell(omega, pop, fluid, x, y);
		}
	}
}

double error(std::size_t x) {
	double acc = 0.0;

	for ( std::size_t y = 1; y < dimY-1; ++y ) {
		acc += std::abs(poiseuille(y) - fluid.velocity(x,y)[0]);
	}

	return acc / (dimY-2);
}

int main() {
	init();

	std::cout << "dim:     " << dimX << "x" << dimY << std::endl;
	std::cout << "Re:      " << reynolds << std::endl;
	std::cout << "uInflow: " << uInflow  << std::endl;
	std::cout << "tau:     " << tau      << std::endl;
	std::cout << "omega:   " << omega    << std::endl;

	std::cout << std::endl;

	for ( std::size_t t = 0; t <= 10000; ++t ) {
		computeLbmStep();

		if ( t % 1000 == 0 ) {
			std::cout << error(dimX-1) << std::endl;
			fluid.writeAsVTK("result/poiseuille_t" + std::to_string(t) + ".vtk");
		}
	}
}
