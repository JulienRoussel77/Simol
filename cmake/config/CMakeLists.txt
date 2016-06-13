# Your project and its supported languages
PROJECT(<your_project_name> CXX)

# The minimum CMake version required
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# Your executable and its dependencies
# Replace <your_executable_name> by the name of your executable
# and my_file.cpp by the list of your source files.
# For example, if you want to generate an executable called foo.exe,
# which depends on main.cpp, bar.cpp and baz.cpp, 
# you should replace the line below by
#
# ADD_EXECUTABLE(foo.exe main.cpp bar.cpp baz.cpp)
#
ADD_EXECUTABLE(<your_executable_name> my_file.cpp)

# Include directories containing header files
INCLUDE_DIRECTORIES(/local/cdoucet/external/yamlcpp/install/include)
INCLUDE_DIRECTORIES(/local/cdoucet/external/eigen/install/include/eigen3)
INCLUDE_DIRECTORIES(/local/cdoucet/my_simol/include)

# Libraries to link with
TARGET_LINK_LIBRARIES(bichain /local/cdoucet/my_simol/lib/libsimol-statphys.a)
TARGET_LINK_LIBRARIES(bichain /local/cdoucet/external/yamlcpp/install/lib/libyaml-cpp.a)

# Compiler flag for C++11 support
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
