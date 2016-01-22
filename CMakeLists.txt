########################
# \file CMakeLists.txt
# \author Cedric Doucet
########################

########################################################################
# Minimum version is 2.8 because it the default version for Ubuntu 12.04
# However, older versions may be used (not tested)
########################################################################

cmake_minimum_required (VERSION 2.8)

################################################################
# The SIMOL_GCC_ROOT has been defined for two main reasons:
#   1. it prevents users from providing paths for each compiler
#   2. it aims at solving linking problems with GCC 5 and higher
# But for the moment, the problem with GCC 5 still remains.
################################################################

if(SIMOL_GCC_ROOT)
  set(CMAKE_C_COMPILER ${SIMOL_GCC_ROOT}/bin/gcc)
  set(CMAKE_CXX_COMPILER ${SIMOL_GCC_ROOT}/bin/g++)
  set(CMAKE_Fortran_COMPILER ${SIMOL_GCC_ROOT}/bin/gfortran)
  link_directories(${SIMOL_GCC_ROOT}/lib64)
endif(SIMOL_GCC_ROOT)

################################################
# Project definition:
#   1. its name is simol
#   2. it involves both fortran, c and c++ files
################################################

project (simol Fortran C CXX)

#==============
# INCLUDE FILES
#==============

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

include(CheckCXXCompilerFlag)
include(ExternalProject)
include(ProcessorCount)
include(SimolInstallPackage)

ProcessorCount(N)
if(NOT N EQUAL 0)
  set(CMAKE_BUILD_FLAGS -j${N})
endif()

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

#=============================
# DOWNLOAD 3RD-PARTY LIBRARIES
#=============================

if(INSTALL_GCC)
  set(GCC_URL ftp://ftp.uvsq.fr/pub/gcc/releases/gcc-4.9.2/gcc-4.9.2.tar.gz)
  set(GCC_CONFIGURE_COMMAND cd <BINARY_DIR> && 
                            <SOURCE_DIR>/contrib/download_prerequisites && 
                            <SOURCE_DIR>/configure --prefix=<INSTALL_DIR>)
#  DownloadThenInstall(GCC)
endif(INSTALL_GCC)


if(ENABLE_DOWNLOAD)
  message(STATUS "libraries to download ${LIBRARIES_TO_DOWNLOAD}")
  simol_install_package(Armadillo)
  simol_install_package(Arpack)
  simol_install_package(Boost)
  #simol_install_package(CppUnit)
  simol_install_package(Eigen)
  simol_install_package(Yamlcpp)
endif(ENABLE_DOWNLOAD)

  include_directories(${EIGEN_INSTALL_DIR}/include/eigen3)
  include_directories(${YAMLCPP_INSTALL_DIR}/include)
  include_directories(${BOOST_INSTALL_DIR}/include)
  include_directories(${ARMADILLO_INSTALL_DIR}/include)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SIMOL_CXX_FLAGS} -Wall -Wextra -fopenmp -lgfortran -Wno-deprecated")

CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
  message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} does not support -std=c++11 flag. Please use a different C++ compiler.")
endif()

if (CMAKE_BUILD_TYPE EQUAL "Release")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2" )
elseif (CMAKE_BUILD_TYPE EQUAL "Profile")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -pg" )
else ()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g" )
endif(CMAKE_BUILD_TYPE EQUAL "Release")

#================
# LINKER SETTINGS
#================


add_subdirectory(src)
