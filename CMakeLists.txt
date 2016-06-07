########################
# Minimum version is 2.8
########################

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

################################################################
# The SIMOL_GCC_ROOT has been defined for two main reasons:
#   1. it prevents users from providing paths for each compiler
#   2. it aims at solving linking problems with GCC 5 and higher
# But for the moment, the problem with GCC 5 still remains.
################################################################

IF(SIMOL_GCC_ROOT)
  SET(CMAKE_C_COMPILER ${SIMOL_GCC_ROOT}/bin/gcc)
  SET(CMAKE_CXX_COMPILER ${SIMOL_GCC_ROOT}/bin/g++)
  SET(CMAKE_Fortran_COMPILER ${SIMOL_GCC_ROOT}/bin/gfortran)
  SET(LD_LIBRARY_PATH ${LD_LIBRARY_PATH} ${SIMOL_GCC_ROOT}/lib64)
  SET(CMAKE_INSTALL_RPATH ${SIMOL_GCC_ROOT}/lib64) 
  LINK_DIRECTORIES(${SIMOL_GCC_ROOT}/lib64)
ENDIF(SIMOL_GCC_ROOT)
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
#SET(CMAKE_EXE_LINKER_FLAGS 
#          "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath -Wl,/usr/local/lib")################################################
# Project definition:
#   1. its name is simol
#   2. it involves both fortran, c and c++ files
################################################

PROJECT(simol Fortran C CXX)

#==================================
# Path to CMake modules
#  - recipes for external libraries
#  - simol functions
#==================================

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

INCLUDE(CheckCXXCompilerFlag)
INCLUDE(ExternalProject)
INCLUDE(SimolInstallPackage)


#===========================
# Automatic setting of make
# to perform parallel builds
#===========================

INCLUDE(ProcessorCount)
ProcessorCount(N)
IF(NOT N EQUAL 0)
  SET(CMAKE_BUILD_FLAGS -j${N})
ENDIF()

#=============
# USER OPTIONS
#=============

OPTION(ENABLE_DOWNLOAD "Enables downloading of missing external libraries" True)
OPTION(INSTALL_GCC "Automatic installation of GCC" False)

IF(NOT EXTERNAL_DIR)
  SET(EXTERNAL_DIR ${CMAKE_SOURCE_DIR}/download)
ENDIF()

FIND_PACKAGE(BLAS)
FIND_PACKAGE(LAPACK)

LIST(GET BLAS_LIBRARIES 0 BLAS_LIB)
LIST(GET LAPACK_LIBRARIES 0 LAPACK_LIB)

#=============================
# DOWNLOAD 3RD-PARTY LIBRARIES
#=============================

IF(INSTALL_GCC)
  SET(GCC_URL ftp://ftp.uvsq.fr/pub/gcc/releases/gcc-4.9.2/gcc-4.9.2.tar.gz)
  SET(GCC_CONFIGURE_COMMAND cd <BINARY_DIR> && 
      <SOURCE_DIR>/contrib/download_prerequisites && 
      <SOURCE_DIR>/configure --prefix=<INSTALL_DIR>)
ENDIF(INSTALL_GCC)


IF(ENABLE_DOWNLOAD)
  MESSAGE(STATUS "libraries to download ${LIBRARIES_TO_DOWNLOAD}")
  simol_install_package(Armadillo)
  simol_install_package(Arpack)
  simol_install_package(Boost)
  simol_install_package(Eigen)
  simol_install_package(Yamlcpp)
ENDIF(ENABLE_DOWNLOAD)

INCLUDE_DIRECTORIES(${EIGEN_INSTALL_DIR}/include/eigen3)
INCLUDE_DIRECTORIES(${YAMLCPP_INSTALL_DIR}/include)
INCLUDE_DIRECTORIES(${BOOST_INSTALL_DIR}/include)
INCLUDE_DIRECTORIES(${ARMADILLO_INSTALL_DIR}/include)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SIMOL_CXX_FLAGS}  -Wall -Wextra -fopenmp -lgfortran -Wno-deprecated")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SIMOL_CXX_FLAGS}  -Werror=return-local-addr -Werror=comment -Werror=sign-compare -Werror=unused-variable -Werror=unused-but-set-variable -Werror=delete-non-virtual-dtor")

CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
IF(COMPILER_SUPPORTS_CXX11)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
ELSE()
  MESSAGE(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} does not support -std=c++11 flag. Please use a different C++ compiler.")
ENDIF()

MESSAGE(STATUS "Build type is ${CMAKE_BUILD_TYPE}")

IF(CMAKE_BUILD_TYPE EQUAL "Release")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3" )
ELSEIF (CMAKE_BUILD_TYPE EQUAL "Profile")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -pg" )
  SET(CMAKE_STATIC_LINKER_FLAGS "${CMAKE_STATIC_LINKER_FLAGS} -pg")
  SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
ELSE()
  # penser a mettre l'option Og pour meilleures perfs en debug
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g" )
ENDIF(CMAKE_BUILD_TYPE EQUAL "Release")

#================
# LINKER SETTINGS
#================

INCLUDE_DIRECTORIES(src)
INCLUDE_DIRECTORIES(include)

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(test)

#=============
# INSTALLATION
#=============
INSTALL(DIRECTORY include/simol 
        DESTINATION include)

