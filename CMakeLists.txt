cmake_minimum_required (VERSION 2.6)

project (simol Fortran C CXX)

#==============
# INCLUDE FILES
#==============

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

include(CheckCXXCompilerFlag)
include(ExternalProject)
include(ProcessorCount)
include(DownloadThenInstall)

set(SIMOL_C_COMPILER ${CMAKE_C_COMPILER})
set(SIMOL_CXX_COMPILER ${CMAKE_CXX_COMPILER})
set(SIMOL_Fortran_COMPILER ${CMAKE_Fortran_COMPILER})

ProcessorCount(N)
if(NOT N EQUAL 0)
  set(CMAKE_BUILD_FLAGS -j${N})
  #set(ctest_test_args ${ctest_test_args} PARALLEL_LEVEL ${N})
endif()
#====================
# EXECUTABLE SETTINGS
#====================
add_executable(test_simol src/test/unit/test_simol.cpp 
               src/test/unit/EigenSolverTest.cpp 
               src/test/unit/SparseMatrixTest.cpp 
               #src/test/unit/MatrixMarketFileTest.cpp 
               src/core/io/mmio.c
               src/test/unit/VectorTest.cpp)

#=============
# USER OPTIONS
#=============

option(ENABLE_DOWNLOAD "Enables downloading of missing external libraries" True)
option(INSTALL_GCC "Automatic installation of GCC" False)

if(NOT EXISTS EXTERNAL_DIR)
  set(EXTERNAL_DIR ${CMAKE_CURRENT_BINARY_DIR}/third-party)
endif()

find_package(BLAS)
find_package(LAPACK)

LIST(GET BLAS_LIBRARIES 0 BLAS_LIB)
LIST(GET LAPACK_LIBRARIES 0 LAPACK_LIB)


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

  include_directories(${EIGEN_INSTALL_DIR}/include/eigen3)
  include_directories(${YAMLCPP_INSTALL_DIR}/include)
  include_directories(${BOOST_INSTALL_DIR}/include)

  message(STATUS "YAMLCPP_INSTALL_DIR is ${YAMLCPP_INSTALL_DIR}")

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

message(STATUS ${CMAKE_LINKER})


#set(SIMOL_CXX_FLAGS "-static")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SIMOL_CXX_FLAGS} -Wall -Wextra -fopenmp -lgfortran -Wno-deprecated")

CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
  message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} does not support -std=c++11 flag. Please use a different C++ compiler.")
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
include_directories(/home/stoltz/travail/Simol/gcc-4.9.2/installation/include)


#set(CMAKE_LINKER /home/stoltz/travail/Simol/gcc-4.9.2/installation/lib64/libstdc++.a)

MESSAGE(STATUS ${BLAS_LIB})
target_link_libraries(test_simol ${CPPUNIT_INSTALL_DIR}/lib/libcppunit.so 
                                 ${LAPACK_LIB}
                                 ${BLAS_LIB}
                                 ${ARPACK_INSTALL_DIR}/lib/libarpack.a
                                 ${ARPACK_INSTALL_DIR}/lib/libarpack_util.a
                                 )

#================
# LINKER SETTINGS
#================

include_directories(src/core)
include_directories(src/core/linalg)
include_directories(src/core/io)
include_directories(src/core/statphys)
include_directories(src/core/random)
include_directories(src)

add_subdirectory(src)

#==============
# TEST SETTINGS
#==============
enable_testing()

add_test(test_simol ${EXECUTABLE_OUTPUT_PATH}/test_simol)

set_tests_properties(test_simol PROPERTIES PASS_REGULAR_EXPRESSION "OK")
