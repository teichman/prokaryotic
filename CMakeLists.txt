cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
project(prokaryotic)
find_package(Eigen3 REQUIRED NO_MODULE)
add_executable(${PROJECT_NAME} prokaryotic.cpp)
target_link_libraries(${PROJECT_NAME} Eigen3::Eigen)


