#pragma once

#include "vector.h"
#include "data_cell_buffer.h"

void computeWallCell(DataCellBuffer& pop, Vector<std::size_t> cell, Vector<int> normal);

void computeMovingWallCell(DataCellBuffer& pop, Vector<std::size_t> cell, Vector<int> normal, Vector<double> u);

void computeZouHeVelocityWallCell(DataCellBuffer& pop, Vector<std::size_t> cell, Vector<int> normal, double vX);
