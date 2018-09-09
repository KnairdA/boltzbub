#pragma once

#include "vector.h"
#include "data_cell_buffer.h"

void computeWallCell(DataCellBuffer& pop, Vector<std::size_t> cell, Vector<int> normal);

void computeZouHeVelocityWallCell(DataCellBuffer& pop, Vector<std::size_t> cell, Vector<int> normal, double vX);
