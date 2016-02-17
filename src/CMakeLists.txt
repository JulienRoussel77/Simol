
# including src/core is useful when one wants to include linalg/Vector.hpp instead of core/linalg/Vector.hpp
#include_directories(core)

ADD_SUBDIRECTORY(core)
ADD_SUBDIRECTORY(quantchem)
ADD_SUBDIRECTORY(statphys)


