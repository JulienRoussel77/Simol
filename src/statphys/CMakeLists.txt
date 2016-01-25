#====================================================
# Store all cpp filenames into cpp-filenames variable
#====================================================
file(GLOB_RECURSE cpp-filenames "*.cpp")

#==============================
# Define simol-statphys library
#==============================
add_library(simol-statphys ${cpp-filenames})

#===================================================
# Force simol-core to be built before simol-statphys
#===================================================
add_dependencies(simol-statphys simol-core)

#===========================================
# Libraries to be linked with simol-statphys
#===========================================
target_link_libraries(simol-statphys simol-core)
target_link_libraries(simol-statphys ${YAMLCPP_INSTALL_DIR}/lib/libyaml-cpp.a)
target_link_libraries(simol-statphys ${ARMADILLO_INSTALL_DIR}/lib/libarmadillo.so)

