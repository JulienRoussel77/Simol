######################
# Compilation settings
######################

SET(SIMOL_COMPILERS_ROOT "/local/cdoucet/gcc/4.9.2" CACHE PATH "Compilers root") # OPTIONAL -- alias for convenience

SET(CMAKE_C_COMPILER "${SIMOL_COMPILERS_ROOT}/gcc" CACHE PATH "C compiler path") # REQUIRED
SET(CMAKE_CXX_COMPILER "${SIMOL_COMPILERS_ROOT}/g++" CACHE PATH "C++ compiler path") # REQUIRED
SET(CMAKE_Fortran_COMPILER "${SIMOL_COMPILERS_ROOT}/gfortran" CACHE PATH "Fortran compiler path") #REQUIRED

#####################
# 3rd party libraries
#####################

SET(SIMOL_AUTO_INSTALL True CACHE BOOL "Automatic downloading and installation of 3rd party libraries") # REQUIRED
SET(SIMOL_EXTERNAL_DIR "/local/cdoucet/external" PATH "External libraries path")

##############
# Installation
##############

SET(SIMOL_INSTALL_PREFIX "/local/cdoucet/my_simol" PATH "Installation path")
