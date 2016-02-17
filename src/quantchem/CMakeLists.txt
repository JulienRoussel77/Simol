#====================================================
# Store all cpp filenames into cpp-filenames variable
#====================================================
FILE(GLOB_RECURSE cpp-filenames "*.cpp")

#==============================
# Define simol-statphys library
#==============================
ADD_LIBRARY(simol-quantchem ${cpp-filenames})

#===================================================
# Force simol-core to be built before simol-statphys
#===================================================
ADD_DEPENDENCIES(simol-quantchem simol-core)

#===========================================
# Libraries to be linked with simol-statphys
#===========================================
TARGET_LINK_LIBRARIES(simol-quantchem simol-core)
TARGET_LINK_LIBRARIES(simol-quantchem ${YAMLCPP_INSTALL_DIR}/lib/libyaml-cpp.a)

#add_executable(quantum_chemistry core/io/mmio.c)
#add_dependencies(quantum_chemistry eigen)
#target_link_libraries(quantum_chemistry ${YAMLCPP_INSTALL_DIR}/lib/libyaml-cpp.a)
#target_link_libraries(quantum_chemistry ${ARPACK_INSTALL_DIR}/lib/libarpack.a)

