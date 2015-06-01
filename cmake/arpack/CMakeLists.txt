
cmake_minimum_required (VERSION 2.6)

project (Arpack Fortran)

#==============
# INCLUDE FILES
#==============

include(ProcessorCount)

ProcessorCount(N)
if(NOT N EQUAL 0)
  set(CMAKE_BUILD_FLAGS -j${N})
  #set(ctest_test_args ${ctest_test_args} PARALLEL_LEVEL ${N})
endif()

#====================
# EXECUTABLE SETTINGS
#====================
file(GLOB SRC "SRC/*.f")

add_library(arpack ${SRC})

