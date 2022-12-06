cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_COMPILE_WARNING_AS_ERROR ON)

project(prokaryotic)

find_package(yaml-cpp REQUIRED)
find_package(Threads REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)
find_package(fmt REQUIRED)
# include(FindProtobuf)
# find_package(Protobuf REQUIRED)
# include_directories(${PROTOBUF_INCLUDE_DIR})
# protobuf_generate_cpp(PROTO_SRC PROTO_HEADER proto/molecule_type.proto)

# macos yaml-cpp cmake entry is broken, obviously.
# fuck this shit though, I'll just build on linux.
# # https://github.com/jbeder/yaml-cpp/issues/774
# # - Config file for the yaml-cpp package
# # It defines the following variables
# #  YAML_CPP_INCLUDE_DIR - include directory
# #  YAML_CPP_LIBRARIES    - libraries to link against
# # Compute paths
# get_filename_component(YAML_CPP_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
# # Our library dependencies (contains definitions for IMPORTED targets)
# include("${YAML_CPP_CMAKE_DIR}/yaml-cpp-targets.cmake")
# # These are IMPORTED targets created by yaml-cpp-targets.cmake
# get_target_property(YAML_CPP_INCLUDE_DIR yaml-cpp INTERFACE_INCLUDE_DIRECTORIES)
# get_target_property(YAML_CPP_LIBRARIES yaml-cpp IMPORTED_LOCATION_RELEASE)

# # This is a really great trick to inspect cmake.  Thanks chatGPT.
get_cmake_property(_variableNames VARIABLES)
foreach(_variableName ${_variableNames})
  message(STATUS "${_variableName}=${${_variableName}}")
endforeach()


include_directories(.)
# include_directories(${PROTOBUF_INCLUDE_DIR})
# include_directories(build/)  # protobuf cc and h files end up here
include_directories(${YAML_CPP_INCLUDE_DIRS})
message("YAML directories: ${YAML_CPP_INCLUDE_DIRS}")

add_compile_options(-fno-omit-frame-pointer -fno-optimize-sibling-calls -O0)
add_link_options(-fno-omit-frame-pointer -fno-optimize-sibling-calls -O0)

# add_compile_options(-fsanitize=address)
# add_link_options(-fsanitize=address)

add_library(prokaryotic SHARED prokaryotic.h prokaryotic.cpp)
target_link_libraries(prokaryotic Eigen3::Eigen fmt::fmt)

add_executable(run prokaryotic.h run.cpp)
target_link_libraries(run prokaryotic Eigen3::Eigen ${YAML_CPP_LIBRARIES})

add_executable(test prokaryotic.h doctest.h test.cpp)
target_link_libraries(test prokaryotic Eigen3::Eigen ${YAML_CPP_LIBRARIES})

# Can't seem to get a stack trace in 
add_executable(debug_snippet prokaryotic.h doctest.h debug_snippet.cpp)
target_link_libraries(debug_snippet prokaryotic Eigen3::Eigen ${YAML_CPP_LIBRARIES})

# add_executable(scratch_protobuf_2 scratch_protobuf_2.cpp ${PROTO_SRC} ${PROTO_HEADER})
# target_link_libraries(scratch_protobuf_2 ${Protobuf_LIBRARIES})

add_executable(scratch_yaml scratch_yaml.cpp)
target_link_libraries(scratch_yaml ${YAML_CPP_LIBRARIES})

