cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_COMPILE_WARNING_AS_ERROR ON)

project(prokaryotic)

find_package(yaml-cpp REQUIRED)
find_package(Threads REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)
find_package(fmt REQUIRED)
find_package(cppzmq)  # https://github.com/zeromq/cppzmq

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
  # This looks like we are mixing library verisons, but:
  #  https://github.com/zeromq/cppzmq
  # "Any libzmq 4.x version is expected to work. DRAFT features may only work for the most recent tested version.
  #  Currently explicitly tested libzmq versions are
  #    4.2.0 (without DRAFT API)
  #    4.3.4 (with and without DRAFT API)"
  set(ZMQ_INCLUDE_DIRS "/opt/homebrew/Cellar/zeromq/4.3.4/include" "/opt/homebrew/Cellar/cppzmq/4.9.0/include/")
  set(ZMQ_LIB_DIRS "/opt/homebrew/Cellar/zeromq/4.3.4/lib")
  set(ZMQ_LIBRARIES "zmq")
endif()

include_directories(${YAML_CPP_INCLUDE_DIRS} ${ZMQ_INCLUDE_DIRS})
link_directories(${YAML_CPP_INCLUDE_DIRS} ${ZMQ_LIB_DIRS})
message("YAML directories: ${YAML_CPP_INCLUDE_DIRS}")

include_directories(.)

add_compile_options(-Wall -Werror)
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  message("Mode: Debug")
  add_compile_options(-fno-omit-frame-pointer -fno-optimize-sibling-calls -O0 -fsanitize=address -g)
  add_link_options(-fno-omit-frame-pointer -fno-optimize-sibling-calls -O0 -fsanitize=address -g)
elseif(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
  message("Mode: RelWithDebInfo")
  add_compile_options(-O3 -g)
  add_link_options(-O3 -g)
endif()

add_library(comms SHARED comms.cpp)
target_link_libraries(comms Eigen3::Eigen fmt::fmt ${YAML_CPP_LIBRARIES} ${ZMQ_LIBRARIES})

add_library(prokaryotic SHARED prokaryotic.h prokaryotic.cpp)
target_link_libraries(prokaryotic Eigen3::Eigen fmt::fmt ${YAML_CPP_LIBRARIES} ${ZMQ_LIBRARIES} comms)

# add_executable(run prokaryotic.h run.cpp)
# target_link_libraries(run prokaryotic Eigen3::Eigen ${YAML_CPP_LIBRARIES})

add_executable(test prokaryotic.h doctest.h test.cpp)
target_link_libraries(test prokaryotic Eigen3::Eigen ${YAML_CPP_LIBRARIES})

add_executable(scratch_yaml scratch_yaml.cpp)
target_link_libraries(scratch_yaml ${YAML_CPP_LIBRARIES})

add_executable(scratch_zmq2python scratch_zmq2python.cpp)
target_link_libraries(scratch_zmq2python comms)

add_executable(scratch_py2cpp scratch_py2cpp.cpp)
target_link_libraries(scratch_py2cpp comms)



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


############################################################
# other crap that is probably done
############################################################

# add_executable(scratch_zmq scratch_zmq.cpp)
# target_link_libraries(scratch_zmq ${ZMQ_LIBRARIES})

# add_executable(debug_snippet prokaryotic.h doctest.h debug_snippet.cpp)
# target_link_libraries(debug_snippet prokaryotic Eigen3::Eigen ${YAML_CPP_LIBRARIES})

