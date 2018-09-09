#pragma once

#include "data_cell_buffer.h"

struct BoxObstacle {
	const std::size_t lower_x_;
	const std::size_t lower_y_;
	const std::size_t upper_x_;
	const std::size_t upper_y_;

	BoxObstacle(std::size_t lX, std::size_t lY, std::size_t uX, std::size_t uY);

	bool isInside(std::size_t x, std::size_t y) const;

	void applyBoundary(DataCellBuffer& pop) const;
};
