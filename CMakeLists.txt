cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_BUILD_TYPE Debug)

find_package(Eigen3 REQUIRED NO_MODULE)

project(prokaryotic)
include_directories(.)

add_library(prokaryotic SHARED prokaryotic.h doctest.h prokaryotic.cpp)
target_link_libraries(prokaryotic Eigen3::Eigen)

add_executable(run prokaryotic.h run.cpp)
target_link_libraries(run prokaryotic Eigen3::Eigen)

add_executable(test prokaryotic.h doctest.h test.cpp)
target_link_libraries(test prokaryotic Eigen3::Eigen)

