#include <iostream>
#include <vector>
#include <algorithm>

#include "lbm.h"
#include "boundary_conditions.h"
#include "box_obstacle.h"

constexpr std::size_t dimX = 500;
constexpr std::size_t dimY = 40;

constexpr double uWall    = 0.2;
constexpr double reynolds = 500;

constexpr double tau   = 3. * uWall * (dimX-1) / reynolds + 0.5;
constexpr double omega = 1. / tau;

DataCellBuffer pop(dimX, dimY);
FluidBuffer fluid(dimX, dimY);

std::vector<BoxObstacle> obstacles{
	{300, 0,  320, 25},
	{340, 15, 360, 39},
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

	// periodic boundary
	for ( std::size_t y = 1; y < dimY-1; ++y ) {
		pop.curr(0,y+1).get(1, 1) = pop.prev(dimX-1,y).get( 1, 1);
		pop.curr(0,y  ).get(1, 0) = pop.prev(dimX-1,y).get( 1, 0);
		pop.curr(0,y-1).get(1,-1) = pop.prev(dimX-1,y).get( 1,-1);
	}
	for ( std::size_t y = 1; y < dimY-1; ++y ) {
		pop.curr(dimX-1,y+1).get(-1, 1) = pop.prev(0,y).get(-1, 1);
		pop.curr(dimX-1,y  ).get(-1, 0) = pop.prev(0,y).get(-1, 0);
		pop.curr(dimX-1,y-1).get(-1,-1) = pop.prev(0,y).get(-1,-1);
	}

	// straight wall cell bounce back
	for ( std::size_t x = 0; x < 100; ++x ) {
		computeZouHeVelocityWallCell(pop, {x, 0     }, { 0, 1}, uWall);
		computeZouHeVelocityWallCell(pop, {x, dimY-1}, { 0,-1}, uWall);
	}

	for ( std::size_t x = 100; x < dimX-1; ++x ) {
		computeWallCell(pop, {x, 0     }, { 0, 1});
		computeWallCell(pop, {x, dimY-1}, { 0,-1});
	}

	// obstacles
	for ( const auto& box : obstacles ) {
		box.applyBoundary(pop);
	}

	for ( std::size_t x = 0; x < dimX; ++x ) {
		for ( std::size_t y = 0; y < dimY; ++y ) {
			Cell& cell = static_cast<Cell&>(pop.curr(x,y));
			fluid.density(x,y)  = cell.sum();
			fluid.velocity(x,y) = cell.velocity(fluid.density(x,y));

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

	std::cout << "Re:    " << reynolds << std::endl;
	std::cout << "uWall: " << uWall    << std::endl;
	std::cout << "tau:   " << tau      << std::endl;
	std::cout << "omega: " << omega    << std::endl;

	for ( std::size_t t = 0; t <= 10000; ++t ) {
		computeLbmStep();

		if ( t % 100 == 0 ) {
			std::cout << ".";
			std::cout.flush();
			fluid.writeAsVTK("result/data_t" + std::to_string(t) + ".vtk");
		}
	}

	std::cout << std::endl;
}
