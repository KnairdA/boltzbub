#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <fstream>

inline double sq(double x) noexcept {
	return x * x;
}

struct DataCell {
	double data[3][3];

	double& get(int x, int y) {
		return data[1+x][1-y];
	}

	double get(int x, int y) const {
		return data[1+x][1-y];
	}
};

struct Vector {
	double data[2];

	double comp(int x, int y) const {
		return x*data[0] + y*data[1];
	}

	double norm() const {
		return std::sqrt(sq(data[0]) + sq(data[1]));
	}

	double& operator[](std::size_t i) {
		return data[i];
	}

	double operator[](std::size_t i) const {
		return data[i];
	}

	Vector operator*(double scalar) const {
		return Vector{
			data[0] * scalar,	
			data[1] * scalar
		};
	}

	Vector& operator+=(const Vector& rhs) {
		data[0] += rhs[0];
		data[1] += rhs[1];
		return *this;
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
		get(-1, 1) = weight.get(-1, 1) * d * ( 1 + 3*v.comp(-1, 1) + 4.5*sq(v.comp(-1, 1)) - 1.5*sq(v.norm()) );
		get( 0, 1) = weight.get( 0, 1) * d * ( 1 + 3*v.comp( 0, 1) + 4.5*sq(v.comp( 0, 1)) - 1.5*sq(v.norm()) );
		get( 1, 1) = weight.get( 1, 1) * d * ( 1 + 3*v.comp( 1, 1) + 4.5*sq(v.comp( 1, 1)) - 1.5*sq(v.norm()) );
		get(-1, 0) = weight.get(-1, 0) * d * ( 1 + 3*v.comp(-1, 0) + 4.5*sq(v.comp(-1, 0)) - 1.5*sq(v.norm()) );
		get( 0, 0) = weight.get( 0, 0) * d * ( 1                                           - 1.5*sq(v.norm()) );
		get( 1, 0) = weight.get( 1, 0) * d * ( 1 + 3*v.comp( 1, 0) + 4.5*sq(v.comp( 1, 0)) - 1.5*sq(v.norm()) );
		get(-1,-1) = weight.get(-1,-1) * d * ( 1 + 3*v.comp(-1,-1) + 4.5*sq(v.comp(-1,-1)) - 1.5*sq(v.norm()) );
		get( 0,-1) = weight.get( 0,-1) * d * ( 1 + 3*v.comp( 0,-1) + 4.5*sq(v.comp( 0,-1)) - 1.5*sq(v.norm()) );
		get( 1,-1) = weight.get( 1,-1) * d * ( 1 + 3*v.comp( 1,-1) + 4.5*sq(v.comp( 1,-1)) - 1.5*sq(v.norm()) );
	}

	double sum() const {
		return get(-1, 1) + get( 0, 1) + get( 1, 1) + get(-1, 0) + get( 0, 0) + get( 1, 0) + get(-1,-1) + get( 0,-1) + get( 1,-1);
	}

	Velocity velocity(Density d) const {
		return Velocity{
			get( 1, 0) - get(-1, 0) + get( 1, 1) - get(-1,-1) + get( 1,-1) - get(-1,1),
			get( 0, 1) - get( 0,-1) + get( 1, 1) - get(-1,-1) - get( 1,-1) + get(-1,1)
		} * ( 1./d );
	}
};

constexpr std::size_t dimX = 128; 
constexpr std::size_t dimY = 128; 

constexpr double tau   = 0.6;
constexpr double omega = 1. / tau;

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
			//velocity[x][y] += force[x][y] * (tau / density[x][y]);

			// compute equilibrium
			Cell eq;
			eq.equilibrize(density[x][y], velocity[x][y]);

			Cell& old = (*oldCells)[x][y];

			// collide & stream
			(*newCells)[x    ][y    ].get( 0, 0) = old.get( 0, 0) + omega * ( eq.get( 0, 0) - old.get( 0, 0) );
			(*newCells)[x + 1][y    ].get( 1, 0) = old.get( 1, 0) + omega * ( eq.get( 1, 0) - old.get( 1, 0) );
			(*newCells)[x - 1][y    ].get(-1, 0) = old.get(-1, 0) + omega * ( eq.get(-1, 0) - old.get(-1, 0) );
			(*newCells)[x    ][y + 1].get( 0, 1) = old.get( 0, 1) + omega * ( eq.get( 0, 1) - old.get( 0, 1) );
			(*newCells)[x    ][y - 1].get( 0,-1) = old.get( 0,-1) + omega * ( eq.get( 0,-1) - old.get( 0,-1) );
			(*newCells)[x + 1][y + 1].get( 1, 1) = old.get( 1, 1) + omega * ( eq.get( 1, 1) - old.get( 1, 1) );
			(*newCells)[x - 1][y - 1].get(-1,-1) = old.get(-1,-1) + omega * ( eq.get(-1,-1) - old.get(-1,-1) );
			(*newCells)[x + 1][y - 1].get( 1,-1) = old.get( 1,-1) + omega * ( eq.get( 1,-1) - old.get( 1,-1) );
			(*newCells)[x - 1][y + 1].get(-1, 1) = old.get(-1, 1) + omega * ( eq.get(-1, 1) - old.get(-1, 1) );
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

	for ( std::size_t x = 50; x < dimX-50; ++x ) {
		for ( std::size_t y = 50; y < dimY-50; ++y ) {
			density[x][y] = 0.8;
		}
	}
}

int main() {
	init();

	for ( std::size_t t = 0; t < 500; ++t ) {
		computeLbmStep(t);

		if ( t % 1 == 0 ) {
			writeCurrentStateAsVTK(t);
		}
	}
}
