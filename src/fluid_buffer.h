#pragma once

#include <memory>
#include <string>
#include <filesystem>

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
	FluidBuffer(std::size_t dimX, std::size_t dimY);

	std::size_t dimX() const;
	std::size_t dimY() const;

	Density& density(std::size_t x, std::size_t y);
	Density& density(Vector<std::size_t> pos);

	Velocity& velocity(std::size_t x, std::size_t y);
	Velocity& velocity(Vector<std::size_t> pos);

	void writeAsVTK(const std::filesystem::path& path);
};
