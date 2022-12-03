cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_COMPILE_WARNING_AS_ERROR ON)

find_package(Eigen3 REQUIRED NO_MODULE)

project(prokaryotic)
include_directories(.)

add_compile_options(-fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls -O0)
add_link_options(-fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls -O0)

add_library(prokaryotic SHARED prokaryotic.h prokaryotic.cpp)
target_link_libraries(prokaryotic Eigen3::Eigen)

add_executable(run prokaryotic.h run.cpp)
target_link_libraries(run prokaryotic Eigen3::Eigen)

add_executable(test prokaryotic.h doctest.h test.cpp)
target_link_libraries(test prokaryotic Eigen3::Eigen)

