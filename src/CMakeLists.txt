
# including src/core is useful when one wants to include linalg/Vector.hpp instead of core/linalg/Vector.hpp
include_directories(core)

add_subdirectory(core)
add_subdirectory(statphys)
add_subdirectory(modules)


