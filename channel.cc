#include <iostream>
#include <vector>
#include <algorithm>

#include "lbm.h"
#include "boundary_conditions.h"
#include "box_obstacle.h"

constexpr std::size_t dimX = 500;
constexpr std::size_t dimY = 40;

constexpr double uInflow  = 0.1;
constexpr double reynolds = 1000;

constexpr double tau   = 3. * uInflow * (dimX-1) / reynolds + 0.5;
constexpr double omega = 1. / tau;

DataCellBuffer pop(dimX, dimY);
FluidBuffer fluid(dimX, dimY);

std::vector<BoxObstacle> obstacles{
	{100, 0,  120, 25},
	{140, 15, 160, 39},
	{300, 15, 340, 25},
};

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
			if ( std::all_of(obstacles.cbegin(), obstacles.cend(), [x, y](const auto& o) {
				return !o.isInside(x, y);
			}) ) {
				streamFluidCell(pop, x, y);
			}
		}
	}

	for ( std::size_t x = 0; x < dimX; ++x ) {
		computeWallCell(pop, {x, 0     }, { 0, 1});
		computeWallCell(pop, {x, dimY-1}, { 0,-1});
	}

	// obstacles
	for ( const auto& box : obstacles ) {
		box.applyBoundary(pop);
	}

	for ( std::size_t x = 0; x < dimX; ++x ) {
		for ( std::size_t y = 0; y < dimY; ++y ) {
			if ( x == 0 && y > 0 && y < dimY-1 ) {
				// inflow
				fluid.density(x,y)  = 1.0;
				fluid.velocity(x,y) = { uInflow, 0 };
			} else {
				Cell& cell = static_cast<Cell&>(pop.curr(x,y));
				fluid.density(x,y)  = cell.sum();
				fluid.velocity(x,y) = cell.velocity(fluid.density(x,y));
			}

			if ( std::all_of(obstacles.cbegin(), obstacles.cend(), [x, y](const auto& o) {
				return !o.isInside(x, y);
			}) ) {
				collideFluidCell(omega, pop, fluid, x, y);
			}
		}
	}
}

int main() {
	init();

	std::cout << "Re:      " << reynolds << std::endl;
	std::cout << "uInflow: " << uInflow  << std::endl;
	std::cout << "tau:     " << tau      << std::endl;
	std::cout << "omega:   " << omega    << std::endl;

	for ( std::size_t t = 0; t <= 5000; ++t ) {
		computeLbmStep();

		if ( t % 50 == 0 ) {
			std::cout << ".";
			std::cout.flush();
			fluid.writeAsVTK("result/data_t" + std::to_string(t) + ".vtk");
		}
	}

	std::cout << std::endl;
}
