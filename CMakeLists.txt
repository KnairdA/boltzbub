cmake_minimum_required(VERSION 3.10)
project(boltzbub LANGUAGES CXX)

add_library(
	boltzbub
		src/lbm.cc
		src/fluid_buffer.cc
		src/boundary_conditions.cc
		src/box_obstacle.cc
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
	lid_driven_cavity
		lid_driven_cavity.cc
)

target_link_libraries(
	lid_driven_cavity
		boltzbub
)

add_executable(
	lid_driven_cavity_with_obstacles
		lid_driven_cavity_with_obstacles.cc
)

target_link_libraries(
	lid_driven_cavity_with_obstacles
		boltzbub
)

add_executable(
	channel
		channel.cc
)

target_link_libraries(
	channel
		boltzbub
)

add_executable(
	poiseuille
		poiseuille.cc
)

target_link_libraries(
	poiseuille
		boltzbub
)
