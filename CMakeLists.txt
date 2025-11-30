##############################################
#              Project Setup                 #
##############################################
# cmake sut up
cmake_minimum_required(VERSION 3.28)
project(01_DD LANGUAGES CXX C)

# c++ standard setup
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O2 ")

# Locate MPI compiler.
find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILER "${MPI_CXX_COMPILER}")

# Locate OpenMP
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP Correctly Initialized")
endif()

# Locate Eigen
find_package (Eigen3 REQUIRED NO_MODULE)

# project structure setup
include_directories(include)
set(01_DD_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${01_DD_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${01_DD_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${01_DD_DIR}/bin")

##############################################
#              Our Libraries                 #
##############################################


##############################################
#                 Our Tests                  #
##############################################
#add_executable()

##############################################
#              Our Executables               #
##############################################
#add_executable()