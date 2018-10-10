cmake_minimum_required(VERSION 3.10)
project(boltzbub)

set(
	CMAKE_CXX_FLAGS
	"${CMAKE_CXX_FLAGS} -std=c++1z -W -Wall -Wextra -Winline -pedantic"
)

include_directories(
	src/
)

set(
	BASE
	src/lbm.cc
	src/fluid_buffer.cc
	src/boundary_conditions.cc
)

add_executable(
	cavity2d
		cavity2d.cc
		${BASE}
)

add_executable(
	cavity2d_with_obstacles
		cavity2d_with_obstacles.cc
		${BASE}
		src/box_obstacle.cc
)
