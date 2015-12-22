add_subdirectory(core)
add_subdirectory(modules)
file(GLOB_RECURSE random "core/random/*.cpp")
file(GLOB_RECURSE io "core/io/*.cpp")
file(GLOB_RECURSE linalg "core/linalg/*.cpp")

add_dependencies(yamlcpp boost)


