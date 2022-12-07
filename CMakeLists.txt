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

# # This is a really great trick to inspect cmake.  Thanks chatGPT.
# get_cmake_property(_variableNames VARIABLES)
# foreach(_variableName ${_variableNames})
#   message(STATUS "${_variableName}=${${_variableName}}")
# endforeach()

# macos yaml-cpp cmake entry is broken, obviously.
if(APPLE)
  # include_directories(/opt/homebrew/Cellar/yaml-cpp/0.7.0/include)
  # link_directories(/opt/homebrew/Cellar/yaml-cpp/0.7.0/lib/)
  set(YAML_CPP_LIBRARIES "yaml-cpp")
  set(YAML_CPP_INCLUDE_DIRS "/opt/homebrew/Cellar/yaml-cpp/0.7.0/include")
endif()

include_directories(${YAML_CPP_INCLUDE_DIRS})
message("YAML directories: ${YAML_CPP_INCLUDE_DIRS}")

include_directories(.)

add_compile_options(-fno-omit-frame-pointer -fno-optimize-sibling-calls -O0)
add_link_options(-fno-omit-frame-pointer -fno-optimize-sibling-calls -O0)

# add_compile_options(-fsanitize=address)
# add_link_options(-fsanitize=address)

add_library(prokaryotic SHARED prokaryotic.h prokaryotic.cpp)
target_link_libraries(prokaryotic Eigen3::Eigen fmt::fmt ${YAML_CPP_LIBRARIES})

# get_target_property(OUT prokaryotic LINK_LIBRARIES)
# message("HELLO ${OUT}")

add_executable(run prokaryotic.h run.cpp)
target_link_libraries(run prokaryotic Eigen3::Eigen ${YAML_CPP_LIBRARIES})

add_executable(test prokaryotic.h doctest.h test.cpp)
target_link_libraries(test prokaryotic Eigen3::Eigen ${YAML_CPP_LIBRARIES})

# add_executable(debug_snippet prokaryotic.h doctest.h debug_snippet.cpp)
# target_link_libraries(debug_snippet prokaryotic Eigen3::Eigen ${YAML_CPP_LIBRARIES})


add_executable(scratch_yaml scratch_yaml.cpp)
target_link_libraries(scratch_yaml ${YAML_CPP_LIBRARIES})



############################################################
# protobuf stuff that I may never return to
############################################################

# add_executable(scratch_protobuf_2 scratch_protobuf_2.cpp ${PROTO_SRC} ${PROTO_HEADER})
# target_link_libraries(scratch_protobuf_2 ${Protobuf_LIBRARIES})

# include_directories(${PROTOBUF_INCLUDE_DIR})
# include_directories(build/)  # protobuf cc and h files end up here

# include(FindProtobuf)
# find_package(Protobuf REQUIRED)
# include_directories(${PROTOBUF_INCLUDE_DIR})
# protobuf_generate_cpp(PROTO_SRC PROTO_HEADER proto/molecule_type.proto)
