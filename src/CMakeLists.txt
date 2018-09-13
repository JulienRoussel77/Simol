#====================================================
# Store all cpp filenames into cpp-filenames variable
#====================================================
FILE(GLOB_RECURSE cpp-filenames *.cpp */*.cpp */*/*.cpp)

INCLUDE_DIRECTORIES(. input output dynamics system potential simulation)

#==============================
# Define simol-statphys library
#==============================
ADD_LIBRARY(simol-statphys ${cpp-filenames})

#===================================================
# Force simol-core to be built before simol-statphys
#===================================================
#ADD_DEPENDENCIES(simol-statphys)

#===========================================
# Libraries to be linked with simol-statphys
#===========================================
TARGET_LINK_LIBRARIES(simol-statphys)
TARGET_LINK_LIBRARIES(simol-statphys ${YAMLCPP_INSTALL_DIR}/lib/libyaml-cpp.a)

INSTALL(TARGETS simol-statphys ARCHIVE DESTINATION lib)

