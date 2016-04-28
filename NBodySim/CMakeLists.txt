CMAKE_MINIMUM_REQUIRED(VERSION 3.1.0 FATAL_ERROR)
cmake_policy(SET CMP0015 NEW)

project(NBodySim CXX)

# set platform defines
if(MSVC)
  add_definitions("-D_CRT_SECURE_NO_WARNINGS" "-DNOMINMAX")
  add_definitions(-D__WINDOWS__)
endif()
if(UNIX)
  add_definitions(-D__LINUX__)
endif()

set (CMAKE_CXX_STANDARD 11)

FIND_PACKAGE(HDF5 NO_MODULE)

include_directories(${HDF5_INCLUDE_DIR})
set (HDF5_libraries hdf5 hdf5_cpp)

add_executable(nbodyrun Simulation.cpp main.cpp)
target_link_libraries(nbodyrun ${HDF5_libraries})

message(STATUS "hello ${HDF5_libraries}")
