
cmake_minimum_required (VERSION 2.6)

project (simol CXX)

#==============
# INCLUDE FILES
#==============

set(CMAKE_VERBOSE_MAKEFILE ON)

include(CheckCXXCompilerFlag)
include(ExternalProject)
include(ProcessorCount)

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


add_subdirectory(src)
#=============
# USER OPTIONS
#=============

option(ENABLE_DOWNLOAD "Enables downloading of missing external libraries" False)
option(EXTERNAL_DIR "Directory where third-party libraries should downloaded" ${CMAKE_CURRENT_BINARY_DIR}/third-party)

#=============================
# DOWNLOAD 3RD-PARTY LIBRARIES
#=============================

if(ENABLE_DOWNLOAD)

  ExternalProject_Add(gcc
                      URL ftp://ftp.uvsq.fr/pub/gcc/releases/gcc-4.9.2/gcc-4.9.2.tar.gz
                      PREFIX ${CMAKE_CURRENT_BINARY_DIR}/third-party/gcc
                      DOWNLOAD_DIR third-party/gcc/download
                      SOURCE_DIR third-party/gcc/configure
                      CONFIGURE_COMMAND cd ../configure && ./contrib/download_prerequisites && ./configure --prefix=${PREFIX}/install
                      #CONFIGURE_COMMAND cd ../gcc && ./contrib/download_prerequisites && ./configure
                      BINARY_DIR third-party/gcc/build
                      BUILD_COMMAND cd ../configure && make
                      INSTALL_DIR third-party/gcc/install
                      INSTALL_COMMAND cd ../configure && make install
  )

#  ExternalProject_Add(yaml-cpp
#                      URL https://yaml-cpp.googlecode.com/files/yaml-cpp-0.5.1.tar.gz
#                      PREFIX ${DOWNLOAD_DIR}/yaml-cpp
#                      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
#  )

ExternalProject_Add(eigen
                    PREFIX ${EXTERNAL_DIR}/eigen
                    DOWNLOAD_DIR ${EXTERNAL_DIR}/eigen/download
                    SOURCE_DIR ${EXTERNAL_DIR}/eigen/src
                    BINARY_DIR ${EXTERNAL_DIR}/eigen/build
                    STAMP_DIR  ${EXTERNAL_DIR}/eigen/stamp
                    INSTALL_DIR ${EXTERNAL_DIR}/eigen/install
                    URL http://bitbucket.org/eigen/eigen/get/3.2.4.tar.gz
                    URL_MD5 4d0d77e06fef87b4fcd2c9b72cc8dc55
                    CMAKE_ARGS -D CMAKE_INSTALL_PREFIX=<INSTALL_DIR> 
                   )

  #ExternalProject_Add(boost
  #                    URL http://sourceforge.net/projects/boost/files/boost/1.57.0/boost_1_57_0.tar.gz
  #                    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/third-party/boost
  #                    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
  #                    CONFIGURE_COMMAND cd ../boost && ./bootstrap.sh
  #                    BUILD_COMMAND cd ../boost && ./b2 install --prefix=.
  #)

  #ExternalProject_Add(cppunit
  #                    URL http://sourceforge.net/projects/cppunit/files/cppunit/1.12.1/cppunit-1.12.1.tar.gz
  #                    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/third-party/cppunit
  #                    CONFIGURE_COMMAND cd ../cppunit && LDFLAGS="-ldl" ./configure --prefix=${PREFIX}
  #)

  add_dependencies(molecular_dynamics gcc)

#  add_dependencies(yaml-cpp boost)

  #add_dependencies(molecular_dynamics cppunit)
#  add_dependencies(molecular_dynamics eigen)
#  add_dependencies(molecular_dynamics yaml-cpp)

#  add_dependencies(quantum_chemistry eigen)
#  add_dependencies(quantum_chemistry yaml-cpp)

  include_directories(${CMAKE_CURRENT_BINARY_DIR}/third-party/eigen/include/eigen3)
  #include_directories(${CMAKE_CURRENT_BINARY_DIR}/third-party/boost/include)
  include_directories(${CMAKE_CURRENT_BINARY_DIR}/third-party/yaml-cpp/include)

  target_link_libraries(molecular_dynamics ${CMAKE_CURRENT_BINARY_DIR}/third-party/yaml-cpp/lib/libyaml-cpp.a)
  target_link_libraries(quantum_chemistry ${CMAKE_CURRENT_BINARY_DIR}/third-party/yaml-cpp/lib/libyaml-cpp.a)

else(ENABLE_DOWNLOAD)

include_directories(/usr/local/bibliotheques/eigen/3.2.4/include/eigen3)
include_directories(/usr/local/bibliotheques/yaml-cpp/0.5.1/include)


target_link_libraries(molecular_dynamics /usr/local/bibliotheques/yaml-cpp/0.5.1/lib/libyaml-cpp.a)
target_link_libraries(quantum_chemistry /usr/local/bibliotheques/yaml-cpp/0.5.1/lib/libyaml-cpp.a)

endif(ENABLE_DOWNLOAD)
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


#set(CMAKE_FORTRAN_COMPILER /home/cdoucet/local/gcc/4.8.2/bin/gfortran)
#set(CMAKE_FORTRAN_COMPILER /usr/local/bibliotheques/g95/g95-install/0.93/bin/x86_64-unknown-linux-gnu-g95)


if (CMAKE_BUILD_TYPE EQUAL "Release")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2" )
elseif (CMAKE_BUILD_TYPE EQUAL "Debug")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g" )
elseif (CMAKE_BUILD_TYPE EQUAL "Profile")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -pg" )
endif(CMAKE_BUILD_TYPE EQUAL "Release")

#==========
# DOWNLOADS
#==========


#====================
# LOCAL INSTALLATIONS
#====================

include_directories(/usr/local/bibliotheques/clapack/3.2.1/include)
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
