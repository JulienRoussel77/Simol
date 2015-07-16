
cmake_minimum_required (VERSION 2.6)

project (simol C CXX)

#==============
# INCLUDE FILES
#==============

#set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

include(CheckCXXCompilerFlag)
include(ExternalProject)
include(ProcessorCount)
include(DownloadThenInstall)

ProcessorCount(N)
if(NOT N EQUAL 0)
  set(CMAKE_BUILD_FLAGS -j${N})
  #set(ctest_test_args ${ctest_test_args} PARALLEL_LEVEL ${N})
endif()
#====================
# EXECUTABLE SETTINGS
#====================
add_executable(quantum_chemistry src/test/functional/quantum_chemistry/main.cpp src/core/io/mmio.c)
add_executable(molecular_dynamics src/test/functional/molecular_dynamics/main.cpp)
add_executable(test_simol src/test/unit/test_simol.cpp src/test/unit/SparseMatrixTest.cpp src/test/unit/MatrixMarketFileTest.cpp src/core/io/mmio.c)


#add_subdirectory(src)
#=============
# USER OPTIONS
#=============

option(ENABLE_DOWNLOAD "Enables downloading of missing external libraries" True)
option(INSTALL_GCC "Automatic installation of GCC" False)

if(NOT EXISTS EXTERNAL_DIR)
  set(EXTERNAL_DIR ${CMAKE_CURRENT_BINARY_DIR}/third-party)
endif()

find_package(Arpack)
find_package(Boost)
find_package(CppUnit)
find_package(Eigen)
find_package(Yamlcpp)
#=============================
# DOWNLOAD 3RD-PARTY LIBRARIES
#=============================

if(INSTALL_GCC)
  set(GCC_URL ftp://ftp.uvsq.fr/pub/gcc/releases/gcc-4.9.2/gcc-4.9.2.tar.gz)
  set(GCC_CONFIGURE_COMMAND cd <BINARY_DIR> && 
                            <SOURCE_DIR>/contrib/download_prerequisites && 
                            <SOURCE_DIR>/configure --prefix=<INSTALL_DIR>)
  DownloadThenInstall(GCC)
endif(INSTALL_GCC)


if(ENABLE_DOWNLOAD)

  message(STATUS "libraries to download ${LIBRARIES_TO_DOWNLOAD}")

  DownloadThenInstall("${LIBRARIES_TO_DOWNLOAD}")



endif(ENABLE_DOWNLOAD)

#  add_dependencies(molecular_dynamics gcc)
  add_dependencies(molecular_dynamics yamlcpp)

  add_dependencies(yamlcpp boost)

  #add_dependencies(molecular_dynamics cppunit)

  add_dependencies(quantum_chemistry eigen)

  include_directories(${EIGEN_INSTALL_DIR}/include/eigen3)
  include_directories(${YAMLCPP_INSTALL_DIR}/include)
  include_directories(${BOOST_INSTALL_DIR}/include)

  message(STATUS "YAMLCPP_INSTALL_DIR is ${YAMLCPP_INSTALL_DIR}")

  target_link_libraries(molecular_dynamics ${YAMLCPP_INSTALL_DIR}/lib/libyaml-cpp.a)
  target_link_libraries(quantum_chemistry ${YAMLCPP_INSTALL_DIR}/lib/libyaml-cpp.a)


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
else()
  message(STATUS "The compiler ${CMAKE_CXX_COMPILER} does not support -std=c++11 flag. Please use a different C++ compiler.")
endif()

if (CMAKE_BUILD_TYPE EQUAL "Release")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2" )
elseif (CMAKE_BUILD_TYPE EQUAL "Debug")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g" )
elseif (CMAKE_BUILD_TYPE EQUAL "Profile")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -pg" )
endif(CMAKE_BUILD_TYPE EQUAL "Release")

#====================
# LOCAL INSTALLATIONS
#====================

include_directories(${CPPUNIT_INSTALL_DIR}/include)
#include_directories(/usr/local/bibliotheques/boost/1.57.0/include)

target_link_libraries(test_simol ${CPPUNIT_INSTALL_DIR}/lib/libcppunit.so ${ARPACK_INSTALL_DIR}/lib/libarpack.a)
target_link_libraries(quantum_chemistry ${ARPACK_INSTALL_DIR}/lib/libarpack.a)

#================
# LINKER SETTINGS
#================

include_directories(src/core)
include_directories(src/core/linalg)
include_directories(src/core/io)
include_directories(src/core/statphys)
include_directories(src)

#==============
# TEST SETTINGS
#==============
enable_testing()

add_test(test_simol ${EXECUTABLE_OUTPUT_PATH}/test_simol)

set_tests_properties(test_simol PROPERTIES PASS_REGULAR_EXPRESSION "OK")
