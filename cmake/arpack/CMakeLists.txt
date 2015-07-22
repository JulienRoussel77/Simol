
cmake_minimum_required (VERSION 2.6)

project (Arpack Fortran)

#==============
# INCLUDE FILES
#==============

include(ProcessorCount)

ProcessorCount(N)
if(NOT N EQUAL 0)
  set(CMAKE_BUILD_FLAGS -j${N})
endif()

#====================
# EXECUTABLE SETTINGS
#====================
file(GLOB SRC "SRC/*.f")
file(GLOB UTIL "UTIL/*.f")

add_library(arpack_util ${UTIL})
add_library(arpack ${SRC})

