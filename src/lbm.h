#pragma once

#include "data_cell.h"
#include "data_cell_buffer.h"
#include "fluid_buffer.h"

constexpr DataCell weight{
	1./36., 1./9., 1./36.,
	1./9.,  4./9., 1./9.,
	1./36., 1./9., 1./36
};

struct Cell : DataCell {
	void equilibrize(Density d, Velocity v);

	double sum() const;

	Velocity velocity(Density d) const;
};

void streamFluidCell(DataCellBuffer& pop, std::size_t x, std::size_t y);

void collideFluidCell(
	double omega,
	DataCellBuffer& pop, FluidBuffer& fluid,
	std::size_t x, std::size_t y);
