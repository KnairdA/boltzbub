#include <iostream>
#include <numeric>
#include <iostream>
#include <fstream>
#include <memory>

#include "src/vector.h"

struct DataCell {
	double data[3][3];

	inline double& get(int x, int y) {
		return data[1+x][1-y];
	}

	inline double& get(Vector<int> v) {
		return get(v[0], v[1]);
	}

	inline double get(int x, int y) const {
		return data[1+x][1-y];
	}

	inline double get(Vector<int> v) const {
		return get(v[0], v[1]);
	}
};

using Velocity = Vector<double>;
using Force = Vector<double>;
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

class CellBuffer {
	private:
		const std::size_t dim_x_;
		const std::size_t dim_y_;

		std::unique_ptr<Cell[]> curr_;
		std::unique_ptr<Cell[]> prev_;

	public:
		CellBuffer(std::size_t dimX, std::size_t dimY):
			dim_x_(dimX),
			dim_y_(dimY),
			curr_(new Cell[dimX*dimY]),
			prev_(new Cell[dimX*dimY]) { }

		void swap() {
			curr_.swap(prev_);
		}

		inline Cell& curr(std::size_t x, std::size_t y) {
			return curr_[y*dim_x_ + x];
		}

		inline Cell& curr(Vector<std::size_t> pos) {
			return curr(pos[0], pos[1]);
		}

		inline Cell& prev(std::size_t x, std::size_t y) {
			return prev_[y*dim_x_ + x];
		}

		inline Cell& prev(Vector<std::size_t> pos) {
			return prev(pos[0], pos[1]);
		}
};

constexpr std::size_t dimX = 128;
constexpr std::size_t dimY = 128;

constexpr double lidVelocityX = 0.05;
constexpr double reynolds     = 1000;
constexpr double nu           = lidVelocityX * dimX / reynolds;
constexpr double tau          = 3*nu + 0.5;

constexpr double omega = 1. / tau;

CellBuffer pop(dimX, dimY);

Density  density [dimX][dimY];
Velocity velocity[dimX][dimY];
Force    force   [dimX][dimY];

void computeFluidCell(std::size_t x, std::size_t y) {
	// compute equilibrium
	Cell eq;
	eq.equilibrize(density[x][y], velocity[x][y]);

	// collide (BGK, relax towards equilibrium) & stream
	for ( int i = -1; i <= 1; ++i ) {
		for ( int j = -1; j <= 1; ++j ) {
			pop.curr(x+i,y+j).get(i,j) = pop.prev(x,y).get(i,j) + omega * (eq.get(i,j) - pop.prev(x,y).get(i,j));
		}
	}
}

std::pair<Vector<int>, Vector<int>> neighbors(Vector<int> v) {
	if ( v[0] == 0 ) {
		return {
			{ -1, v[1] },
			{  1, v[1] }
		};
	} else if ( v[1] == 0 ) {
		return {
			{ v[0], -1 },
			{ v[0],  1 }
		};
	} else {
		return {
			{    0, v[1] },
			{ v[0], 0    }
		};
	}
}

void computeWallCell(Vector<std::size_t> cell, Vector<int> normal) {
	const auto [neighborA, neighborB] = neighbors(normal);

	pop.curr(cell).get(neighborA) = pop.curr(cell-neighborA).get(-neighborA);
	pop.curr(cell).get(normal   ) = pop.curr(cell-normal   ).get(-normal   );
	pop.curr(cell).get(neighborB) = pop.curr(cell-neighborB).get(-neighborB);
}

void computeSimpleVelocityWallCell(Vector<std::size_t> cell, Vector<int> normal, double vX) {
	const auto [neighborA, neighborB] = neighbors(normal);

	pop.curr(cell).get(neighborA) = pop.curr(cell-neighborA).get(-neighborA) - 6 * weight.get( 1, 1) * density[cell[0]][cell[1]] * vX;
	pop.curr(cell).get(normal   ) = pop.curr(cell-normal   ).get(-normal   );
	pop.curr(cell).get(neighborB) = pop.curr(cell-neighborB).get(-neighborB) + 6 * weight.get(-1, 1) * density[cell[0]][cell[0]] * vX;
}

void computeZouHeVelocityWallCell(Vector<std::size_t> cell, Vector<int> normal, double vX) {
	const auto [neighborA, neighborB] = neighbors(normal);

	const double rho = pop.curr(cell-neighborA).get(-1,0) + pop.curr(cell-normal).get(0,0) + pop.curr(cell-neighborB).get(1,0)
		+ 2.*(
			pop.curr(cell-neighborA).get(-neighborA) +
			pop.curr(cell-normal   ).get(-normal   ) +
			pop.curr(cell-neighborB).get(-neighborB)
		);
	const double hor = pop.curr(cell-neighborA).get(-1,0) - pop.curr(cell-neighborB).get(1,0);

	pop.curr(cell).get(neighborA) = pop.curr(cell-neighborA).get(-neighborA) + 0.5*hor - 0.5*vX*rho;
	pop.curr(cell).get(normal   ) = pop.curr(cell-normal   ).get(-normal   );
	pop.curr(cell).get(neighborB) = pop.curr(cell-neighborB).get(-neighborB) - 0.5*hor + 0.5*vX*rho;
}

void computeLbmStep(std::size_t t) {
	pop.swap();

	for ( std::size_t x = 1; x < dimX - 1; ++x ) {
		for ( std::size_t y = 1; y < dimY - 1; ++y ) {
			//if ( x <= 20 || x >= 40 || y <= 20 || y >= 80 ) {
				computeFluidCell(x, y);
			//}
		}
	}

	// obstacle
	/*{
		for ( std::size_t x = 21; x < 40; ++x ) {
			computeWallCell(x, 80, { 0, 1});
			computeWallCell(x, 20, { 0,-1});
		}
		for ( std::size_t y = 21; y < 80; ++y ) {
			computeWallCell(40, y, { 1, 0});
			computeWallCell(20, y, {-1, 0});
		}

		computeWallCell(40, 80, { 1, 1});
		computeWallCell(40, 20, { 1,-1});
		computeWallCell(20, 80, {-1, 1});
		computeWallCell(20, 20, {-1,-1});
	}*/

	// straight wall cell bounce back
	for ( std::size_t x = 2; x < dimX - 2; ++x ) {
		computeWallCell({x, 1}, { 0, 1});
		computeZouHeVelocityWallCell({x, dimY-2}, { 0,-1}, lidVelocityX);
	}
	for ( std::size_t y = 2; y < dimY - 2; ++y ) {
		computeWallCell({1,      y}, { 1, 0});
		computeWallCell({dimX-2, y}, {-1, 0});
	}

	// edge wall cell bounce back
	computeWallCell({1,      1     }, { 1, 1});
	computeWallCell({1,      dimY-2}, { 1,-1});
	computeWallCell({dimX-2, 1     }, {-1, 1});
	computeWallCell({dimX-2, dimY-2}, {-1,-1});

	// update density, velocity field
	for ( std::size_t x = 0; x < dimX; ++x ) {
		for ( std::size_t y = 0; y < dimY; ++y ) {
			density[x][y]  = pop.curr(x,y).sum();
			velocity[x][y] = pop.curr(x,y).velocity(density[x][y]);
		}
	}
}

void init() {
	for ( std::size_t x = 0; x < dimX; ++x ) {
		for ( std::size_t y = 0; y < dimY; ++y ) {
			density[x][y] = 1.0;
			velocity[x][y] = { 0.0, 0.0 };
			force[x][y] = { 0.0, 0.0 };

			pop.curr(x,y).equilibrize(density[x][y], velocity[x][y]);
			pop.prev(x,y).equilibrize(density[x][y], velocity[x][y]);
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

int main() {
	init();

	std::cout << "u: " << lidVelocityX << std::endl;
	std::cout << "tau: " << tau << std::endl;

	for ( std::size_t t = 0; t < 20000; ++t ) {
		computeLbmStep(t);

		if ( t % 1000 == 0 ) {
			std::cout << "Writing " << t << "." << std::endl;
			writeCurrentStateAsVTK(t);
		}
	}
}
