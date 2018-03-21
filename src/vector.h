#include <cmath>

inline double sq(double x) noexcept {
	return x * x;
}

struct Vector {
	double data[2];

	double comp(int x, int y) const {
		return x*data[0] + y*data[1];
	}

	double norm() const {
		return std::sqrt(sq(data[0]) + sq(data[1]));
	}

	double& operator[](std::size_t i) {
		return data[i];
	}

	double operator[](std::size_t i) const {
		return data[i];
	}

	Vector operator*(double scalar) const {
		return Vector{
			data[0] * scalar,
			data[1] * scalar
		};
	}

	Vector& operator+=(const Vector& rhs) {
		data[0] += rhs[0];
		data[1] += rhs[1];
		return *this;
	}
};

Vector operator*(double scalar, const Vector& vec) {
	return Vector{
		vec[0] * scalar,
		vec[1] * scalar
	};
}
