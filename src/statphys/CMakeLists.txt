#====================================================
# Store all cpp filenames into cpp-filenames variable
#====================================================
FILE(GLOB_RECURSE cpp-filenames *.cpp */*.cpp */*/*.cpp)

INCLUDE_DIRECTORIES(. input output controlVariate dynamics system potential system/chain dynamics/boundaryLangevin dynamics/dpde)

#==============================
# Define simol-statphys library
#==============================
ADD_LIBRARY(simol-statphys ${cpp-filenames})

#===================================================
# Force simol-core to be built before simol-statphys
#===================================================
ADD_DEPENDENCIES(simol-statphys simol-core)

#===========================================
# Libraries to be linked with simol-statphys
#===========================================
TARGET_LINK_LIBRARIES(simol-statphys simol-core)
TARGET_LINK_LIBRARIES(simol-statphys ${YAMLCPP_INSTALL_DIR}/lib/libyaml-cpp.a)
TARGET_LINK_LIBRARIES(simol-statphys ${ARMADILLO_INSTALL_DIR}/lib/libarmadillo.so)

