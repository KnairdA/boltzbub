cmake_minimum_required(VERSION 3.10)
project(boltzbub)

set(
	CMAKE_CXX_FLAGS
	"${CMAKE_CXX_FLAGS} -std=c++1z -W -Wall -Wextra -Winline -pedantic"
)

include_directories(
	src/
)

add_executable(
	boltz
		boltz.cc
		src/lbm.cc
		src/fluid_buffer.cc
		src/boundary_conditions.cc
		src/box_obstacle.cc
)
