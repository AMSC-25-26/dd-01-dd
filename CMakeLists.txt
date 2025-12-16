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
find_package(OpenMP REQUIRED)
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
add_library(tridiag src/tridiag_utils.cpp)

add_library(solver src/solver.cpp)
target_link_libraries(solver PRIVATE OpenMP::OpenMP_CXX)

add_custom_target(library)
add_dependencies(library tridiag solver)

##############################################
#                 Our Tests                  #
##############################################
add_executable(test_tridiag test/test_tridiag.cpp)
target_link_libraries(test_tridiag tridiag OpenMP::OpenMP_CXX)

add_executable(test_solver test/test_solver.cpp)
target_link_libraries(test_solver solver tridiag OpenMP::OpenMP_CXX)

add_executable(test_solver_highfreq test/test_solver_highfreq.cpp)
target_link_libraries(test_solver_highfreq solver tridiag OpenMP::OpenMP_CXX)

add_executable(test_solver_reaction test/test_solver_reaction.cpp)
target_link_libraries(test_solver_reaction solver tridiag OpenMP::OpenMP_CXX)

add_executable(test_solver_polynomial test/test_solver_polynomial.cpp)
target_link_libraries(test_solver_polynomial solver tridiag OpenMP::OpenMP_CXX)

add_executable(test_solver_shifted_domain test/test_solver_shifted_domain.cpp)
target_link_libraries(test_solver_shifted_domain solver tridiag OpenMP::OpenMP_CXX)

add_executable(test_solver_overlap_steep test/test_solver_overlap_steep.cpp)
target_link_libraries(test_solver_overlap_steep solver tridiag OpenMP::OpenMP_CXX)

add_custom_target(tests)
add_dependencies(tests
        library
        test_tridiag
        test_solver
        test_solver_highfreq
        test_solver_reaction
        test_solver_polynomial
        test_solver_shifted_domain
        test_solver_overlap_steep
)


##############################################
#              Our Executables               #
##############################################
#add_executable()


##############################################
#              Our Documentation             #
##############################################
find_package(Doxygen)
if(Doxygen_FOUND)
    message(STATUS "Doxygen found")
    doxygen_add_docs(docs ${01_DD_DIR} WORKING_DIRECTORY ${01_DD_DIR} CONFIG_FILE "${01_DD_DIR}/docs/Doxygen/Doxyfile")
endif()