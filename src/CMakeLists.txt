add_subdirectory(core)

file(GLOB_RECURSE src src/*.cpp)

add_executable(quantum_chemistry test/functional/quantum_chemistry/main.cpp core/io/mmio.c)
add_executable(molecular_dynamics test/functional/molecular_dynamics/main.cpp ${src})

add_dependencies(quantum_chemistry eigen)
target_link_libraries(quantum_chemistry ${YAMLCPP_INSTALL_DIR}/lib/libyaml-cpp.a)
target_link_libraries(quantum_chemistry ${ARPACK_INSTALL_DIR}/lib/libarpack.a)


add_dependencies(molecular_dynamics yamlcpp)
add_dependencies(yamlcpp boost)
target_link_libraries(molecular_dynamics ${YAMLCPP_INSTALL_DIR}/lib/libyaml-cpp.a)

