#pragma once

#include <memory>

#include "data_cell.h"

class DataCellBuffer {
private:
	const std::size_t dim_x_;
	const std::size_t dim_y_;

	std::unique_ptr<DataCell[]> curr_;
	std::unique_ptr<DataCell[]> prev_;

public:
	DataCellBuffer(std::size_t dimX, std::size_t dimY):
		dim_x_(dimX),
		dim_y_(dimY),
		curr_(new DataCell[dimX*dimY]),
		prev_(new DataCell[dimX*dimY]) { }

	std::size_t dimX() const {
		return dim_x_;
	}

	std::size_t dimY() const {
		return dim_y_;
	}

	void swap() {
		curr_.swap(prev_);
	}

	inline DataCell& curr(std::size_t x, std::size_t y) {
		return curr_[y*dim_x_ + x];
	}

	inline DataCell& curr(Vector<std::size_t> pos) {
		return curr(pos[0], pos[1]);
	}

	inline DataCell& prev(std::size_t x, std::size_t y) {
		return prev_[y*dim_x_ + x];
	}

	inline DataCell& prev(Vector<std::size_t> pos) {
		return prev(pos[0], pos[1]);
	}
};

