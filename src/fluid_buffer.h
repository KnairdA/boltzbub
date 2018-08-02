#pragma once

#include <fstream>
#include <memory>

#include "vector.h"

using Velocity = Vector<double>;
using Density  = double;

class FluidBuffer {
	private:
		const std::size_t dim_x_;
		const std::size_t dim_y_;

		std::unique_ptr<Density[]> density_;
		std::unique_ptr<Velocity[]> velocity_;

	public:
		FluidBuffer(std::size_t dimX, std::size_t dimY):
			dim_x_(dimX),
			dim_y_(dimY),
			density_(new Density[dimX*dimY]),
			velocity_(new Velocity[dimX*dimY]) { }

		std::size_t dimX() const {
			return dim_x_;
		}

		std::size_t dimY() const {
			return dim_y_;
		}

		Density& density(std::size_t x, std::size_t y) {
			return density_[y*dim_x_ + x];
		}

		Density& curr(Vector<std::size_t> pos) {
			return density(pos[0], pos[1]);
		}

		Velocity& velocity(std::size_t x, std::size_t y) {
			return velocity_[y*dim_x_ + x];
		}

		Velocity& velocity(Vector<std::size_t> pos) {
			return velocity(pos[0], pos[1]);
		}

		void writeAsVTK(int time) {
			std::ofstream fout;
			fout.open(("result/data_t" + std::to_string(time) + ".vtk").c_str());

			fout << "# vtk DataFile Version 3.0\n";
			fout << "lbm_output\n";
			fout << "ASCII\n";
			fout << "DATASET RECTILINEAR_GRID\n";
			fout << "DIMENSIONS " << dimX() << " " << dimY() << " 1" << "\n";

			fout << "X_COORDINATES " << dimX() << " float\n";
			for( std::size_t x = 0; x < dimX(); ++x ) {
				fout << x << " ";
			}

			fout << "\nY_COORDINATES " << dimY() << " float\n";
			for( std::size_t y = 0; y < dimY(); ++y ) {
				fout << y << " ";
			}

			fout << "\nZ_COORDINATES " << 1 << " float\n";
			fout << 0 << "\n";
			fout << "POINT_DATA " << dimX() * dimY() << "\n";

			fout << "VECTORS velocity float\n";
			for ( std::size_t y = 0; y < dimY(); ++y ) {
				for ( std::size_t x = 0; x < dimX(); ++x ) {
					fout << velocity(x,y)[0] << " " << velocity(x,y)[1] << " 0\n";
				}
			}

			fout << "SCALARS density float 1\n";
			fout << "LOOKUP_TABLE default\n";
			for ( std::size_t y = 0; y < dimY(); ++y ) {
				for ( std::size_t x = 0; x < dimX(); ++x ) {
					fout << density(x,y) << "\n";
				}
			}

			fout.close();
		}
};
