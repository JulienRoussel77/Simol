
# including src/core is useful when one wants to include linalg/Vector.hpp instead of core/linalg/Vector.hpp
include_directories(core)
include_directories(core/linalg)
include_directories(core/io)
include_directories(core/statphys)
include_directories(core/random)

add_subdirectory(core)
add_subdirectory(modules)


