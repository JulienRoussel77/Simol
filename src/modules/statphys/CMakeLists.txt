file(GLOB_RECURSE statphys "modules/statphys/*.cpp")

#add_executable(molecular_dynamics ${statphys} ${random} ${io} ${linalg})



#add_dependencies(molecular_dynamics armadillo)
#add_dependencies(molecular_dynamics yamlcpp)

#target_link_libraries(molecular_dynamics ${ARMADILLO_INSTALL_DIR}/lib/libarmadillo.so)
#target_link_libraries(molecular_dynamics ${YAMLCPP_INSTALL_DIR}/lib/libyaml-cpp.a)

