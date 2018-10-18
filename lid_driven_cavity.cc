#include <iostream>

#include "lbm.h"
#include "boundary_conditions.h"

constexpr std::size_t dimX = 100;
constexpr std::size_t dimY = dimX;

constexpr double uLid     = 0.4;
constexpr double reynolds = 1000;

constexpr double tau   = 3. * uLid * (dimX-1) / reynolds + 0.5;
constexpr double omega = 1. / tau;

DataCellBuffer pop(dimX, dimY);
FluidBuffer fluid(dimX, dimY);

void init() {
	for ( std::size_t x = 0; x < dimX; ++x ) {
		for ( std::size_t y = 0; y < dimY; ++y ) {
			fluid.density(x,y)  = 1.0;
			fluid.velocity(x,y) = { 0.0, 0.0 };

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

	// moving top wall
	for ( std::size_t x = 0; x < dimX; ++x ) {
		computeMovingWallCell(pop, {x, dimY-1}, {0, -1}, {uLid, 0});
	}

	// straight wall cell bounce back
	for ( std::size_t x = 1; x < dimX-1; ++x ) {
		computeWallCell(pop, {x, 0}, { 0, 1});
	}
	for ( std::size_t y = 1; y < dimY-1; ++y ) {
		computeWallCell(pop, {0,      y}, { 1, 0});
		computeWallCell(pop, {dimX-1, y}, {-1, 0});
	}

	// edge wall cell bounce back
	computeWallCell(pop, {0,      0 }, { 1, 1});
	computeWallCell(pop, {dimX-1, 0 }, {-1, 1});

	for ( std::size_t x = 0; x < dimX; ++x ) {
		for ( std::size_t y = 0; y < dimY; ++y ) {
			Cell& cell = static_cast<Cell&>(pop.curr(x,y));
			fluid.density(x,y)  = cell.sum();
			fluid.velocity(x,y) = cell.velocity(fluid.density(x,y));

			collideFluidCell(omega, pop, fluid, x, y);
		}
	}
}

int main() {
	init();

	std::cout << "Re:    " << reynolds << std::endl;
	std::cout << "uLid:  " << uLid     << std::endl;
	std::cout << "tau:   " << tau      << std::endl;
	std::cout << "omega: " << omega    << std::endl;

	for ( std::size_t t = 0; t <= 10000; ++t ) {
		computeLbmStep();

		if ( t % 1000 == 0 ) {
			std::cout << ".";
			std::cout.flush();
			fluid.writeAsVTK("result/lid_driven_cavity_t" + std::to_string(t) + ".vtk");
		}
	}

	std::cout << std::endl;
}
