cmake_minimum_required(VERSION 3.10)
project(boltzbub LANGUAGES CXX)

add_library(
	boltzbub
		src/lbm.cc
		src/fluid_buffer.cc
		src/boundary_conditions.cc
)

target_compile_features(
	boltzbub
	PUBLIC
		cxx_std_17
)

target_include_directories(
	boltzbub
	PUBLIC
		${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_link_libraries(
	boltzbub
	PRIVATE
		stdc++fs
)

add_executable(
	cavity2d
		cavity2d.cc
)

target_link_libraries(
	cavity2d
		boltzbub
)

add_executable(
	cavity2d_with_obstacles
		cavity2d_with_obstacles.cc
		src/box_obstacle.cc
)

target_link_libraries(
	cavity2d_with_obstacles
		boltzbub
)
