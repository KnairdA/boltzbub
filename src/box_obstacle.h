#pragma once

#include "data_cell_buffer.h"
#include "boundary_conditions.h"

struct BoxObstacle {
	const std::size_t lower_x_;
	const std::size_t lower_y_;
	const std::size_t upper_x_;
	const std::size_t upper_y_;

	BoxObstacle(std::size_t lX, std::size_t lY, std::size_t uX, std::size_t uY):
		lower_x_(lX), lower_y_(lY), upper_x_(uX), upper_y_(uY) { }

	bool isInside(std::size_t x, std::size_t y) const {
		return x > lower_x_
		    && x < upper_x_
		    && y > lower_y_
		    && y < upper_y_;
	}

	void applyBoundary(DataCellBuffer& pop) const {
		for ( std::size_t x = lower_x_+1; x < upper_x_; ++x ) {
			computeWallCell(pop, {x, lower_y_}, { 0,-1});
			computeWallCell(pop, {x, upper_y_}, { 0, 1});
		}
		for ( std::size_t y = lower_y_+1; y < upper_y_; ++y ) {
			computeWallCell(pop, {lower_x_, y}, {-1, 0});
			computeWallCell(pop, {upper_x_, y}, { 1, 0});
		}
		computeWallCell(pop, {lower_x_, lower_y_}, {-1,-1});
		computeWallCell(pop, {upper_x_, lower_y_}, { 1,-1});
		computeWallCell(pop, {upper_x_, upper_y_}, { 1, 1});
		computeWallCell(pop, {lower_x_, upper_y_}, {-1, 1});
	}
};
