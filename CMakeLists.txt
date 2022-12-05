cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_COMPILE_WARNING_AS_ERROR ON)

project(prokaryotic)

find_package(Threads REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)
include(FindProtobuf)
find_package(Protobuf REQUIRED)
include_directories(${PROTOBUF_INCLUDE_DIR})
protobuf_generate_cpp(PROTO_SRC PROTO_HEADER proto/addressbook.proto)

include_directories(.)
include_directories(${PROTOBUF_INCLUDE_DIR})
include_directories(build/)  # protobuf cc and h files end up here

add_compile_options(-fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls -O0)
add_link_options(-fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls -O0)

add_library(prokaryotic SHARED prokaryotic.h prokaryotic.cpp)
target_link_libraries(prokaryotic Eigen3::Eigen)

add_executable(run prokaryotic.h run.cpp)
target_link_libraries(run prokaryotic Eigen3::Eigen)

add_executable(test prokaryotic.h doctest.h test.cpp ${PROTO_SRC} ${PROTO_HEADER})
target_link_libraries(test prokaryotic Eigen3::Eigen ${Protobuf_LIBRARIES})

add_executable(scratch_protobuf scratch_protobuf.cpp ${PROTO_SRC} ${PROTO_HEADER})
target_link_libraries(scratch_protobuf ${Protobuf_LIBRARIES})


