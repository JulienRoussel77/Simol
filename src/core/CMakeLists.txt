file(GLOB_RECURSE random "random/*.cpp")
file(GLOB_RECURSE io "io/*.cpp")
file(GLOB_RECURSE linalg "linalg/*.cpp")

add_dependencies(yamlcpp boost)

add_library(simol-core SHARED ${random} ${io} ${linalg})

add_custom_target(core DEPENDS simol-core)
add_dependencies(core eigen)
include_directories(${EIGEN_INSTALL_DIR}/include/eigen3)

