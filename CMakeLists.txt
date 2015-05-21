
cmake_minimum_required (VERSION 2.6)

project (simol CXX)

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
add_executable(quantum_chemistry src/test/functional/quantum_chemistry/main.cpp)
add_executable(molecular_dynamics src/test/functional/molecular_dynamics/main.cpp)
add_executable(test_simol src/test/unit/test_simol.cpp)


#add_subdirectory(src)
#=============
# USER OPTIONS
#=============

option(ENABLE_DOWNLOAD "Enables downloading of missing external libraries" False)

if(NOT EXISTS EXTERNAL_DIR)
  set(EXTERNAL_DIR ${CMAKE_CURRENT_BINARY_DIR}/third-party)
endif()

find_package(Boost)
find_package(CppUnit)
find_package(Eigen)
find_package(Yamlcpp)
#=============================
# DOWNLOAD 3RD-PARTY LIBRARIES
#=============================

if(ENABLE_DOWNLOAD)

  message(STATUS "libraries to download ${LIBRARIES_TO_DOWNLOAD}")

  DownloadThenInstall("${LIBRARIES_TO_DOWNLOAD}")

#  ExternalProject_Add(gcc
#                      URL ftp://ftp.uvsq.fr/pub/gcc/releases/gcc-4.9.2/gcc-4.9.2.tar.gz
#                      PREFIX ${CMAKE_CURRENT_BINARY_DIR}/third-party/gcc
#                      DOWNLOAD_DIR third-party/gcc/download
#                      SOURCE_DIR third-party/gcc/configure
#                      CONFIGURE_COMMAND cd ../configure && ./contrib/download_prerequisites && ./configure --prefix=${PREFIX}/install
                      #CONFIGURE_COMMAND cd ../gcc && ./contrib/download_prerequisites && ./configure
#                      BINARY_DIR third-party/gcc/build
#                      BUILD_COMMAND cd ../configure && make
#                      INSTALL_DIR third-party/gcc/install
#                      INSTALL_COMMAND cd ../configure && make install
#  )


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

target_link_libraries(test_simol ${CPPUNIT_INSTALL_DIR}/lib/libcppunit.so)
target_link_libraries(quantum_chemistry /usr/local/bibliotheques/arpack/2.1/lib/libarpack_LINUX.a)

#================
# LINKER SETTINGS
#================

include_directories(src/core)
include_directories(src/core/statphys)
include_directories(src)

#==============
# TEST SETTINGS
#==============
enable_testing()

add_test(test_simol ${EXECUTABLE_OUTPUT_PATH}/test_simol)

set_tests_properties(test_simol PROPERTIES PASS_REGULAR_EXPRESSION "OK")
