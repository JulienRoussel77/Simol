#include_directories(linalg)
#include_directories(io)
include_directories(statphys)
#include_directories(random)

file(GLOB_RECURSE random "random/*.cpp")
file(GLOB_RECURSE io "io/*.cpp")
file(GLOB_RECURSE linalg "linalg/*.cpp")

add_dependencies(yamlcpp boost)

add_library(simol-core SHARED ${random} ${io} ${linalg})

add_custom_target(core DEPENDS simol-core)
add_dependencies(simol-core eigen)

