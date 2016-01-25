file(GLOB_RECURSE statphys "*.cpp")

#add_executable(molecular_dynamics ${statphys} ${random} ${io} ${linalg})
#add_executable(simol-statphys ${statphys})
add_library(simol-statphys ${statphys})
add_dependencies(simol-statphys simol-core)

message(STATUS ${CMAKE_BINARY_DIR})

#target_link_libraries(simol-statphys ${CMAKE_BINARY_DIR}/src/core/simol-core.so)
target_link_libraries(simol-statphys simol-core)

#add_dependencies(molecular_dynamics armadillo)
add_dependencies(simol-statphys yamlcpp)
target_link_libraries(simol-statphys ${YAMLCPP_INSTALL_DIR}/lib/libyaml-cpp.a)
#target_link_libraries(simol-statphys yamlcpp)

target_link_libraries(simol-statphys ${ARMADILLO_INSTALL_DIR}/lib/libarmadillo.so)
#target_link_libraries(molecular_dynamics ${YAMLCPP_INSTALL_DIR}/lib/libyaml-cpp.a)

