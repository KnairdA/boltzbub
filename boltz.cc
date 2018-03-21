#include <iostream>
#include <numeric>
#include <iostream>
#include <fstream>

#include "src/vector.h"

struct DataCell {
	double data[3][3];

	inline double& get(int x, int y) {
		return data[1+x][1-y];
	}

	inline double get(int x, int y) const {
		return data[1+x][1-y];
	}
};

using Velocity = Vector;
using Force = Vector;
using Density = double;

constexpr DataCell weight{
	1./36., 1./9., 1./36.,
	1./9.,  4./9., 1./9.,
	1./36., 1./9., 1./36
};

struct Cell : DataCell {
	void equilibrize(Density d, Velocity v) {
		for ( int i = -1; i <= 1; ++i ) {
			for ( int j = -1; j <= 1; ++j ) {
				get(i,j) = weight.get(i,j) * d * (1 + 3*v.comp(i,j) + 4.5*sq(v.comp(i,j)) - 1.5*sq(v.norm()));
			}
		}
	}

	double sum() const {
		return get(-1, 1) + get( 0, 1) + get( 1, 1) + get(-1, 0) + get( 0, 0) + get( 1, 0) + get(-1,-1) + get( 0,-1) + get( 1,-1);
	}

	Velocity velocity(Density d) const {
		return 1./d * Velocity{
			get( 1, 0) - get(-1, 0) + get( 1, 1) - get(-1,-1) + get( 1,-1) - get(-1,1),
			get( 0, 1) - get( 0,-1) + get( 1, 1) - get(-1,-1) - get( 1,-1) + get(-1,1)
		};
	}
};

constexpr std::size_t dimX = 128;
constexpr std::size_t dimY = 128;

constexpr double tau = 1.0;

Cell cellsA[dimX][dimY];
Cell cellsB[dimX][dimY];

Cell (*newCells)[dimX][dimY] = &cellsA;
Cell (*oldCells)[dimX][dimY] = &cellsB;

Density  density [dimX][dimY];
Velocity velocity[dimX][dimY];
Force    force   [dimX][dimY];

void swapCellBuffers() {
	if ( newCells == &cellsA && oldCells == &cellsB ) {
		newCells = &cellsB;
		oldCells = &cellsA;
	} else {
		newCells = &cellsA;
		oldCells = &cellsB;
	}
}

void computeLbmStep(std::size_t t) {
	swapCellBuffers();

	for ( std::size_t x = 1; x < dimX - 1; ++x ) {
		for ( std::size_t y = 1; y < dimY - 1; ++y ) {
			// update velocity (force coupling)
			//velocity[x][y] += tau / density[x][y] * force[x][y];

			// compute equilibrium
			Cell eq;
			eq.equilibrize(density[x][y], velocity[x][y]);

			// collide & stream
			for ( int i = -1; i <= 1; ++i ) {
				for ( int j = -1; j <= 1; ++j ) {
					(*newCells)[x+i][y+j].get(i,j) = (*oldCells)[x][y].get(i,j) + 1./tau * (eq.get(i,j) - (*oldCells)[x][y].get(i,j));
				}
			}
		}
	}

	// straight wall cell bounce back
	for ( std::size_t x = 2; x < dimX - 2; ++x ) {
		(*newCells)[x][1].get( 0, 1) = (*newCells)[x][0].get( 0,-1);
		(*newCells)[x][1].get( 1, 1) = (*newCells)[x][0].get(-1,-1);
		(*newCells)[x][1].get(-1, 1) = (*newCells)[x][0].get( 1,-1);

		(*newCells)[x][dimY-2].get( 0,-1) = (*newCells)[x][dimY-2].get( 0, 1);
		(*newCells)[x][dimY-2].get(-1,-1) = (*newCells)[x][dimY-2].get( 1, 1);
		(*newCells)[x][dimY-2].get( 1,-1) = (*newCells)[x][dimY-2].get(-1, 1);
	}
	for ( std::size_t y = 2; y < dimY - 2; ++y ) {
		(*newCells)[1][y].get( 1, 0) = (*newCells)[0][y].get(-1, 0);
		(*newCells)[1][y].get( 1,-1) = (*newCells)[0][y].get(-1, 1);
		(*newCells)[1][y].get( 1, 1) = (*newCells)[0][y].get(-1,-1);

		(*newCells)[dimX-2][y].get(-1, 0) = (*newCells)[dimX - 1][y].get( 1, 0);
		(*newCells)[dimX-2][y].get(-1, 1) = (*newCells)[dimX - 1][y].get( 1,-1);
		(*newCells)[dimX-2][y].get(-1,-1) = (*newCells)[dimX - 1][y].get( 1, 1);
	}

	// edge wall cell bounce back
	(*newCells)[1][1].get( 0, 1) = (*newCells)[0][0].get( 0,-1);
	(*newCells)[1][1].get( 1, 1) = (*newCells)[0][0].get(-1,-1);
	(*newCells)[1][1].get( 1, 0) = (*newCells)[0][0].get(-1, 0);

	(*newCells)[1][dimY-2].get( 0,-1) = (*newCells)[0][dimY-1].get( 0, 1);
	(*newCells)[1][dimY-2].get( 1,-1) = (*newCells)[0][dimY-1].get(-1, 1);
	(*newCells)[1][dimY-2].get( 1, 0) = (*newCells)[0][dimY-1].get(-1, 0);

	(*newCells)[dimX-2][1].get( 0, 1) = (*newCells)[dimX-1][0].get( 0,-1);
	(*newCells)[dimX-2][1].get(-1, 1) = (*newCells)[dimX-1][0].get( 1,-1);
	(*newCells)[dimX-2][1].get(-1, 0) = (*newCells)[dimX-1][0].get( 1, 0);

	(*newCells)[dimX-2][dimY-2].get( 0,-1) = (*newCells)[dimX-1][dimY-1].get( 0, 1);
	(*newCells)[dimX-2][dimY-2].get(-1,-1) = (*newCells)[dimX-1][dimY-1].get( 1, 1);
	(*newCells)[dimX-2][dimY-2].get(-1, 0) = (*newCells)[dimX-1][dimY-1].get( 1, 0);

	// update density, velocity field
	for ( std::size_t x = 0; x < dimX; ++x ) {
		for ( std::size_t y = 0; y < dimY; ++y ) {
			density[x][y]  = (*newCells)[x][y].sum();
			velocity[x][y] = (*newCells)[x][y].velocity(density[x][y]);
		}
	}
}

void writeCurrentStateAsVTK(int time) {
	std::ofstream fout;
	fout.open(("result/data_t" + std::to_string(time) + ".vtk").c_str());

	fout << "# vtk DataFile Version 3.0\n";
	fout << "lbm_output\n";
	fout << "ASCII\n";
	fout << "DATASET RECTILINEAR_GRID\n";
	fout << "DIMENSIONS " << dimX - 2 << " " << dimY - 2 << " 1" << "\n";

	fout << "X_COORDINATES " << dimX - 2 << " float\n";
	for( std::size_t x = 1; x < dimX - 1; ++x ) {
		fout << x << " ";
	}

	fout << "\nY_COORDINATES " << dimY - 2 << " float\n";
	for( std::size_t y = 1; y < dimY - 1; ++y ) {
		fout << y << " ";
	}

	fout << "\nZ_COORDINATES " << 1 << " float\n";
	fout << 0 << "\n";
	fout << "POINT_DATA " << (dimX - 2) * (dimY - 2) << "\n";

	fout << "VECTORS velocity float\n";
	for ( std::size_t y = 1; y < dimY - 1; ++y ) {
		for ( std::size_t x = 1; x < dimX - 1; ++x ) {
			fout << velocity[x][y][0] << " " << velocity[x][y][1] << " 0\n";
		}
	}

	fout << "SCALARS density float 1\n";
	fout << "LOOKUP_TABLE default\n";
	for ( std::size_t y = 1; y < dimY - 1; ++y ) {
		for ( std::size_t x = 1; x < dimX - 1; ++x ) {
			fout << density[x][y] << "\n";
		}
	}

	fout.close();
}

void init() {
	for ( std::size_t x = 0; x < dimX; ++x ) {
		for ( std::size_t y = 0; y < dimY; ++y ) {
			density[x][y] = 1.0;
			velocity[x][y] = { 0.0, 0.0 };
			force[x][y] = { 0.0, 0.0 };

			(*newCells)[x][y].equilibrize(density[x][y], velocity[x][y]);
			(*oldCells)[x][y].equilibrize(density[x][y], velocity[x][y]);
		}
	}

	for ( std::size_t y = 55; y < dimY-55; ++y ) {
		for ( std::size_t x = 35; x < dimX-75; ++x ) {
			density[x][y] = 0.8;
		}
		for ( std::size_t x = 75; x < dimX-35; ++x ) {
			density[x][y] = 0.8;
		}
	}
}

int main() {
	init();

	for ( std::size_t t = 0; t < 100; ++t ) {
		computeLbmStep(t);

		if ( t % 2 == 0 ) {
			writeCurrentStateAsVTK(t);
		}
	}
}
