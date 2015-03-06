
cmake_minimum_required (VERSION 2.6)

project (simol CXX Fortran)

#==============
# INCLUDE FILES
#==============

include(CheckCXXCompilerFlag)
include(ExternalProject)

#====================
# EXECUTABLE SETTINGS
#====================
add_executable(quantum_chemistry test/functional/quantum_chemistry/main.cpp)
add_executable(molecular_dynamics test/functional/molecular_dynamics/main.cpp)

#=============================
# DOWNLOAD 3RD-PARTY LIBRARIES
#=============================

ExternalProject_Add(yaml-cpp
                    URL https://yaml-cpp.googlecode.com/files/yaml-cpp-0.5.1.tar.gz
                    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/third-party/yaml-cpp
                    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
)

add_dependencies(molecular_dynamics yaml-cpp)

#===================
# USER CONFIGURATION
#===================

set (config ${CMAKE_CURRENT_SOURCE_DIR}/config.cmake)
if(EXISTS ${config})
  include(${config})
endif(EXISTS ${config})

#==================
# COMPILER SETTINGS
#==================

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -fopenmp")

CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else(COMPILER_SUPPORTS_CXX11)
      #  message(STATUS "The compiler ${CMAKE_CXX_COMPILER} does not support -std=c++11 flag. Please use a different C++ compiler.")
endif(COMPILER_SUPPORTS_CXX11)


#set(CMAKE_FORTRAN_COMPILER /home/cdoucet/local/gcc/4.8.2/bin/gfortran)
#set(CMAKE_FORTRAN_COMPILER /usr/local/bibliotheques/g95/g95-install/0.93/bin/x86_64-unknown-linux-gnu-g95)


if (CMAKE_BUILD_TYPE EQUAL "Release")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2" )
elseif (CMAKE_BUILD_TYPE EQUAL "Debug")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g" )
elseif (CMAKE_BUILD_TYPE EQUAL "Profile")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -pg" )
endif(CMAKE_BUILD_TYPE EQUAL "Release")

#====================
# INCLUDE DIRECTORIES
#====================

include_directories(/usr/local/bibliotheques/eigen/3.2.3/include/eigen3)
include_directories(/usr/local/bibliotheques/clapack/3.2.1/include)
include_directories(/usr/local/bibliotheques/cppunit/1.12.1/include)
#include_directories(/usr/local/bibliotheques/yaml-cpp/0.5.1/include)
include_directories(/usr/local/bibliotheques/boost/1.57.0/include)

include_directories(${CMAKE_CURRENT_BINARY_DIR}/third-party/yaml-cpp/include)
include_directories(src/statistical_physics)
include_directories(src)


#================
# LINKER SETTINGS
#================

target_link_libraries(quantum_chemistry /usr/local/bibliotheques/yaml-cpp/0.5.1/lib/libyaml-cpp.a)
target_link_libraries(quantum_chemistry /usr/local/bibliotheques/arpack/2.1/lib/libarpack_LINUX.a)
target_link_libraries(quantum_chemistry /usr/local/bibliotheques/clapack/3.2.1/lib/lapack_LINUX.a)
target_link_libraries(quantum_chemistry cblas)

target_link_libraries(molecular_dynamics ${CMAKE_CURRENT_BINARY_DIR}/third-party/yaml-cpp/lib/libyaml-cpp.a)
#target_link_libraries(molecular_dynamics /usr/local/bibliotheques/yaml-cpp/0.5.1/lib/libyaml-cpp.a)

#==============
# TEST SETTINGS
#==============
enable_testing()
add_executable(test_simol test/unit/test_simol.cpp)
target_link_libraries(test_simol /usr/local/bibliotheques/cppunit/1.12.1/lib/libcppunit.so)

add_test(test_simol ${EXECUTABLE_OUTPUT_PATH}/test_simol)

set_tests_properties(test_simol PROPERTIES PASS_REGULAR_EXPRESSION "OK")
