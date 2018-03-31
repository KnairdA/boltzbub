#include <cmath>

inline double sq(double x) noexcept {
	return x * x;
}

template <typename T>
struct Vector {
	T data[2];

	T comp(int x, int y) const {
		return x*data[0] + y*data[1];
	}

	T norm() const {
		return std::sqrt(sq(data[0]) + sq(data[1]));
	}

	T& operator[](std::size_t i) {
		return data[i];
	}

	T operator[](std::size_t i) const {
		return data[i];
	}

	Vector<T> operator*(T scalar) const {
		return Vector<T>{
			data[0] * scalar,
			data[1] * scalar
		};
	}

	Vector<T>& operator+=(const Vector<T>& rhs) {
		data[0] += rhs[0];
		data[1] += rhs[1];
		return *this;
	}
};

template <typename T>
Vector<T> operator*(T scalar, const Vector<T>& vec) {
	return Vector<T>{
		vec[0] * scalar,
		vec[1] * scalar
	};
}
