#include <iostream>
#include <numeric>
#include <iostream>
#include <fstream>
#include <memory>

#include "vector.h"

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

constexpr std::size_t dimX = 256;
constexpr std::size_t dimY = dimX;

constexpr double uLid     = 0.1;
constexpr double reynolds = 1000;

constexpr double tau   = 3 * uLid * (dimX-1) / reynolds + 0.5;
constexpr double omega = 1. / tau;

CellBuffer pop(dimX, dimY);

Density  density [dimX][dimY];
Velocity velocity[dimX][dimY];

void streamFluidCell(std::size_t x, std::size_t y) {
	if ( x != 0 && x != dimX - 1 && y != 0 && y != dimY -1 ) {
		// stream internal cells
		for ( int i = -1; i <= 1; ++i ) {
			for ( int j = -1; j <= 1; ++j ) {
				pop.curr(x+i,y+j).get(i,j) = pop.prev(x,y).get(i,j);
			}
		}
	} else {
		// stream boundary cells,
		// missing populations to be determined by boundary conditions
		for ( int i = -1; i <= 1; ++i ) {
			for ( int j = -1; j <= 1; ++j ) {
				if ( x+i >= 0 && x+i <= dimX - 1 && y+j >= 0 && y+j <= dimY - 1 ) {
					pop.curr(x+i,y+j).get(i,j) = pop.prev(x,y).get(i,j);
				}
			}
		}
	}
}

void collideFluidCell(std::size_t x, std::size_t y) {
	// compute equilibrium
	Cell eq;
	eq.equilibrize(density[x][y], velocity[x][y]);

	// collide (BGK, relax towards equilibrium)
	if ( x != 0 && x != dimX - 1 && y != 0 && y != dimY -1 ) {
		for ( int i = -1; i <= 1; ++i ) {
			for ( int j = -1; j <= 1; ++j ) {
				pop.curr(x,y).get(i,j) = pop.curr(x,y).get(i,j) + omega * (eq.get(i,j) - pop.curr(x,y).get(i,j));
			}
		}
	} else {
		// partial collide for boundary cells
		for ( int i = -1; i <= 1; ++i ) {
			for ( int j = -1; j <= 1; ++j ) {
				if ( x+i >= 0 && x+i <= dimX - 1 && y+j >= 0 && y+j <= dimY - 1 ) {
					pop.curr(x,y).get(i,j) = pop.curr(x,y).get(i,j) + omega * (eq.get(i,j) - pop.curr(x,y).get(i,j));
				}
			}
		}
	}
}

void computeWallCell(Vector<std::size_t> cell, Vector<int> normal) {
	const auto [neighborA, neighborB] = neighbors(normal);

	pop.curr(cell).get(neighborA) = pop.curr(cell).get(-neighborA);
	pop.curr(cell).get(normal   ) = pop.curr(cell).get(-normal   );
	pop.curr(cell).get(neighborB) = pop.curr(cell).get(-neighborB);
}

void computeZouHeVelocityWallCell(Vector<std::size_t> cell, Vector<int> normal, double vX) {
	const auto [neighborA, neighborB] = neighbors(normal);

	const double rho = pop.curr(cell).get(-1,0) + pop.curr(cell).get(0,0) + pop.curr(cell).get(1,0)
		+ 2.*(
			pop.curr(cell).get(-neighborA) +
			pop.curr(cell).get(-normal   ) +
			pop.curr(cell).get(-neighborB)
		);

	pop.curr(cell).get(neighborA) = pop.curr(cell).get(-neighborA) + 0.5*(pop.curr(cell).get( 1,0) - pop.curr(cell).get(-1,0) - vX*rho);
	pop.curr(cell).get(normal   ) = pop.curr(cell).get(-normal   );
	pop.curr(cell).get(neighborB) = pop.curr(cell).get(-neighborB) + 0.5*(pop.curr(cell).get(-1,0) - pop.curr(cell).get( 1,0) + vX*rho);
}

void computeLbmStep(std::size_t t) {
	pop.swap();

	for ( std::size_t x = 0; x < dimX; ++x ) {
		for ( std::size_t y = 0; y < dimY; ++y ) {
			streamFluidCell(x, y);
		}
	}

	// straight wall cell bounce back
	for ( std::size_t x = 0; x < dimX; ++x ) {
		computeZouHeVelocityWallCell({x, dimY-1}, { 0,-1}, uLid);
	}
	for ( std::size_t x = 1; x < dimX-1; ++x ) {
		computeWallCell({x, 0}, { 0, 1});
	}
	for ( std::size_t y = 1; y < dimY-1; ++y ) {
		computeWallCell({0,      y}, { 1, 0});
		computeWallCell({dimX-1, y}, {-1, 0});
	}

	// edge wall cell bounce back
	computeWallCell({0,      0 }, { 1, 1});
	computeWallCell({dimX-1, 0 }, {-1, 1});

	for ( std::size_t x = 0; x < dimX; ++x ) {
		for ( std::size_t y = 0; y < dimY; ++y ) {
			density[x][y]  = pop.curr(x,y).sum();
			velocity[x][y] = pop.curr(x,y).velocity(density[x][y]);

			collideFluidCell(x, y);
		}
	}
}

void init() {
	for ( std::size_t x = 0; x < dimX; ++x ) {
		for ( std::size_t y = 0; y < dimY; ++y ) {
			velocity[x][y] = { 0.0, 0.0 };
			density[x][y] = 1.0;

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
	fout << "DIMENSIONS " << dimX << " " << dimY << " 1" << "\n";

	fout << "X_COORDINATES " << dimX << " float\n";
	for( std::size_t x = 0; x < dimX; ++x ) {
		fout << x << " ";
	}

	fout << "\nY_COORDINATES " << dimY << " float\n";
	for( std::size_t y = 0; y < dimY; ++y ) {
		fout << y << " ";
	}

	fout << "\nZ_COORDINATES " << 1 << " float\n";
	fout << 0 << "\n";
	fout << "POINT_DATA " << dimX * dimY << "\n";

	fout << "VECTORS velocity float\n";
	for ( std::size_t y = 0; y < dimY; ++y ) {
		for ( std::size_t x = 0; x < dimX; ++x ) {
			fout << velocity[x][y][0] << " " << velocity[x][y][1] << " 0\n";
		}
	}

	fout << "SCALARS density float 1\n";
	fout << "LOOKUP_TABLE default\n";
	for ( std::size_t y = 0; y < dimY; ++y ) {
		for ( std::size_t x = 0; x < dimX; ++x ) {
			fout << density[x][y] << "\n";
		}
	}

	fout.close();
}

int main() {
	init();

	std::cout << "Re:   " << reynolds << std::endl;
	std::cout << "uLid: " << uLid     << std::endl;
	std::cout << "tau:  " << tau      << std::endl;

	for ( std::size_t t = 0; t <= 100000; ++t ) {
		computeLbmStep(t);

		if ( (t < 1000 && t % 10 == 0) || t % 1000 == 0 ) {
			std::cout << ".";
			std::cout.flush();
			writeCurrentStateAsVTK(t);
		}
	}

	std::cout << std::endl;
}
