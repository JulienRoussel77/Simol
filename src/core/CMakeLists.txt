INCLUDE_DIRECTORIES(.)


FILE(GLOB_RECURSE random "random/*.cpp")
FILE(GLOB_RECURSE io "io/*.cpp" "io/*.c")
FILE(GLOB_RECURSE linalg "linalg/*.cpp")

ADD_DEPENDENCIES(yamlcpp boost)

ADD_LIBRARY(simol-core SHARED ${random} ${io} ${linalg})

ADD_CUSTOM_TARGET(core DEPENDS simol-core)
ADD_DEPENDENCIES(simol-core eigen)

