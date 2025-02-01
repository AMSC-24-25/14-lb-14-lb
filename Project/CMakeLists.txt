cmake_minimum_required(VERSION 3.10)

# Set project name
project(MyProject)

# Explicitly use mpic++ as the compiler
find_package(MPI REQUIRED)
if(MPI_CXX_FOUND)
    message(STATUS "MPI found: ${MPI_CXX_COMPILER}")
    set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})  # Force mpic++ as compiler
    include_directories(${MPI_INCLUDE_PATH})
endif()

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set compiler flags (matching your manual command)
set(CMAKE_CXX_FLAGS "-std=c++17 -pedantic -O3 -Wall -fopenmp")

# Add include directories
include_directories(/usr/include ${CMAKE_SOURCE_DIR})

# List source files explicitly
set(SRC_FILES
    LBM.cpp
    VelocitySet.cpp
    LidDrivenCavity3D.cpp
    ObstacleLiftDrag.cpp
    main.cpp
)

# Find and enable OpenMP
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found!")
endif()

# Create the executable
add_executable(ldc_sim_mpi ${SRC_FILES})

# Link OpenMP and MPI correctly
target_link_libraries(ldc_sim_mpi PUBLIC OpenMP::OpenMP_CXX MPI::MPI_CXX)
