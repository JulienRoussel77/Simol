######################
# Compilation settings
######################

#SET(SIMOL_COMPILERS_ROOT "path where looking for compilers" CACHE PATH "Compilers root") # OPTIONAL -- alias for convenience

#SET(CMAKE_C_COMPILER "${SIMOL_COMPILERS_ROOT}/bin/gcc" CACHE PATH "C compiler path") # REQUIRED
#SET(CMAKE_CXX_COMPILER "${SIMOL_COMPILERS_ROOT}/bin/g++" CACHE PATH "C++ compiler path") # REQUIRED
#SET(CMAKE_Fortran_COMPILER "${SIMOL_COMPILERS_ROOT}/bin/gfortran" CACHE PATH "Fortran compiler path") #REQUIRED

#####################
# 3rd party libraries
#####################

SET(SIMOL_AUTO_INSTALL True CACHE BOOL "Automatic downloading and installation of 3rd party libraries") # REQUIRED
SET(SIMOL_EXTERNAL_DIR "/Users/ci/external" CACHE PATH "External libraries path")

##############
# Installation
##############

SET(SIMOL_INSTALL_PREFIX "/Users/ci/installed_simol" CACHE PATH "Installation path")
