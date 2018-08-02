#pragma once

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
