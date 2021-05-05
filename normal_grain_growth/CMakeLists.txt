cmake_minimum_required(VERSION 3.5.1)

project(grain_growth)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE Release)
cmake_policy(SET CMP0071 NEW)

# get the root folder of digitalTwin
get_filename_component(PROJECT_SOURCE_DIR ${CMAKE_SOURCE_DIR} DIRECTORY)
message(STATUS "Project root folder: ${PROJECT_SOURCE_DIR}")

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)

FIND_PACKAGE(MPI REQUIRED)
MESSAGE("\nFound MPI: ${MPI_USE_FILE}")
include_directories(${MPI_INCLUDE_PATH})

find_package(VTK 8.2 REQUIRED HINTS /opt/vtk-8.2)
if (VTK_FOUND)
	message(STATUS "${BoldGreen}Found VTK: ${VTK_USE_FILE}${ColourReset}")
	include(${VTK_USE_FILE})
else()
    message(FATAL_ERROR "${BoldRed}VTK not found${ColourReset}")
endif()	

file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/screenshots)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(HEADER_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${HEADER_DIR})

include_directories(".")
subdirs(src)
