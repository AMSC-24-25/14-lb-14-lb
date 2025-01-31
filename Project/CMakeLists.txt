cmake_minimum_required(VERSION 3.10)
project(LatticeBoltzmannMethod)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3")

file(GLOB SRC_FILES src/*.cpp)
add_library(lbm ${SRC_FILES})
target_include_directories(lbm
                           PUBLIC
                                "${CMAKE_SOURCE_DIR}/include"
                                "/usr/include/eigen3")

find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
        message(STATUS "OpenMP found!")
        target_compile_options(lbm PRIVATE ${OpenMP_CXX_FLAGS})
        target_link_options(lbm PUBLIC ${OpenMP_CXX_FLAGS})
elseif(NOT OpenMP_CXX_FOUND)
        error("OpenMp not found")
endif()

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
if(Eigen3_FOUND)
    message("Eigen Found")
    target_include_directories(lbm PRIVATE ${Eigen_INCLUDE_DIRS})
    target_link_libraries(lbm PRIVATE Eigen3::Eigen)
elseif(NOT Eigen3_FOUND)
    error("Eigen Not Found")
endif()

find_package(MPI REQUIRED)
if(MPI_FOUND)
        message("MPI Found")
        target_link_libraries(lbm PUBLIC MPI::MPI_CXX)
elseif(NOT MPI_FOUND)
        error("MPI Not Found")
endif()

get_property(include_dirs TARGET lbm PROPERTY INCLUDE_DIRECTORIES)
message(STATUS "Include Directories: ${include_dirs}")

add_executable(lbm_sim main.cpp)
target_link_libraries(lbm_sim lbm)