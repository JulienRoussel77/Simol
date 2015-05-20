
cmake_minimum_required (VERSION 2.6)

project (simol CXX)

#==============
# INCLUDE FILES
#==============

set(CMAKE_VERBOSE_MAKEFILE ON)

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


find_package(Eigen)
#find_package(Yamlcpp)
#=============================
# DOWNLOAD 3RD-PARTY LIBRARIES
#=============================

if(ENABLE_DOWNLOAD)

  DownloadThenInstall(${LIBRARIES_TO_DOWNLOAD})

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

#  ExternalProject_Add(yaml-cpp
#                      PREFIX ${EXTERNAL_DIR}/yaml-cpp
#                      DOWNLOAD_DIR ${EXTERNAL_DIR}/yaml-cpp/download
#                      SOURCE_DIR ${EXTERNAL_DIR}/yaml-cpp/src
#                      BINARY_DIR ${EXTERNAL_DIR}/yaml-cpp/build
#                      STAMP_DIR  ${EXTERNAL_DIR}/yaml-cpp/stamp
#                      INSTALL_DIR ${EXTERNAL_DIR}/yaml-cpp/install
#                      URL https://yaml-cpp.googlecode.com/files/yaml-cpp-0.5.1.tar.gz
#                      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
#  )


#  ExternalProject_Add(boost
#                      PREFIX ${EXTERNAL_DIR}/boost
#                      DOWNLOAD_DIR ${EXTERNAL_DIR}/boost/download
#                      SOURCE_DIR ${EXTERNAL_DIR}/boost/src
#                      BINARY_DIR ${EXTERNAL_DIR}/boost/src
#                      STAMP_DIR  ${EXTERNAL_DIR}/boost/stamp
#                      INSTALL_DIR ${EXTERNAL_DIR}/boost/install
#                      URL http://sourceforge.net/projects/boost/files/boost/1.57.0/boost_1_57_0.tar.gz
#                      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
#                      CONFIGURE_COMMAND ./bootstrap.sh
#                      BUILD_COMMAND ./b2 install --prefix=<INSTALL_DIR>
#  )

#  ExternalProject_Add(cppunit
#                      PREFIX ${EXTERNAL_DIR}/third-party/cppunit
#                      SOURCE_DIR ${EXTERNAL_DIR}/cppunit/src
#                      BINARY_DIR ${EXTERNAL_DIR}/cppunit/build
#                      STAMP_DIR  ${EXTERNAL_DIR}/cppunit/stamp
#                      INSTALL_DIR ${EXTERNAL_DIR}/cppunit/install
#                      URL http://sourceforge.net/projects/cppunit/files/cppunit/1.12.1/cppunit-1.12.1.tar.gz
#                      CONFIGURE_COMMAND cd ../cppunit && LDFLAGS="-ldl" ./configure --prefix=${PREFIX}
#  )

endif(ENABLE_DOWNLOAD)

  add_dependencies(molecular_dynamics gcc)
  add_dependencies(molecular_dynamics yamlcpp)

#  add_dependencies(yaml-cpp boost)

  #add_dependencies(molecular_dynamics cppunit)
#  add_dependencies(molecular_dynamics eigen)

#  add_dependencies(quantum_chemistry eigen)
#  add_dependencies(quantum_chemistry yaml-cpp)

  message(STATUS "EIGEN_INSTALL_DIR is ${EIGEN_INSTALL_DIR}")

  include_directories(${EIGEN_INSTALL_DIR}/include/eigen3)
  #include_directories(${CMAKE_CURRENT_BINARY_DIR}/third-party/boost/include)
  include_directories(${YAMLCPP_INSTALL_DIR}/include)

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

include_directories(/usr/local/bibliotheques/cppunit/1.12.1/include)
include_directories(/usr/local/bibliotheques/boost/1.57.0/include)

target_link_libraries(test_simol /usr/local/bibliotheques/cppunit/1.12.1/lib/libcppunit.so)
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
