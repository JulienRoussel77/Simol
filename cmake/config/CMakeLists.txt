# Your project and its supported languages
PROJECT(main-project CXX)

# The minimum CMake version required
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# Your executable and its dependencies
# Replace <your-executable-name> by the name of your executable
# and my_file.cpp by the list of your source files.
# For example, if you want to generate an executable called foo.exe,
# which depends on main.cpp, bar.cpp and baz.cpp, 
# you should replace the line below by
#
# ADD_EXECUTABLE(foo.exe main.cpp bar.cpp baz.cpp)
#
ADD_EXECUTABLE(main main.cpp)

# Include directories containing header files
INCLUDE_DIRECTORIES(${SIMOL_EXTERNAL_DIR}/yamlcpp/install/include)
INCLUDE_DIRECTORIES(${SIMOL_EXTERNAL_DIR}/eigen/install/include/eigen3)
INCLUDE_DIRECTORIES(${SIMOL_INSTALL_PREFIX}/include)

# Libraries to link with
TARGET_LINK_LIBRARIES(main ${SIMOL_INSTALL_PREFIX}/lib/libsimol-core.so)
TARGET_LINK_LIBRARIES(main ${SIMOL_INSTALL_PREFIX}/lib/libsimol-quantchem.a)
TARGET_LINK_LIBRARIES(main ${SIMOL_INSTALL_PREFIX}/lib/libsimol-statphys.a)
TARGET_LINK_LIBRARIES(main ${SIMOL_EXTERNAL_DIR}/yamlcpp/install/lib/libyaml-cpp.a)

# Compiler flag for C++11 support
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")
