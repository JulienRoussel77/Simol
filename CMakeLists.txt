
cmake_minimum_required (VERSION 2.6)

project (simol)

#====================
# EXECUTABLE SETTINGS
#====================
add_executable(simol src/main.cpp)

#==================
# COMPILER SETTINGS
#==================

set(CMAKE_CXX_COMPILER /home/cdoucet/local/gcc/4.8.2/bin/g++)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -fopenmp -std=c++11")


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

#================
# LINKER SETTINGS
#================

#link_directories(/usr/local/bibliotheques/clapack/3.2.1/lib)
target_link_libraries(simol cblas)
target_link_libraries(simol /usr/local/bibliotheques/clapack/3.2.1/lib/lapack_LINUX.a)


#==============
# TEST SETTINGS
#==============
enable_testing()
add_executable(test_simol test/unit/test_simol.cpp)
target_link_libraries(test_simol /usr/local/bibliotheques/cppunit/1.12.1/lib/libcppunit.so)

add_test(test_simol ${EXECUTABLE_OUTPUT_PATH}/test_simol)

set_tests_properties(test_simol PROPERTIES PASS_REGULAR_EXPRESSION "OK")
