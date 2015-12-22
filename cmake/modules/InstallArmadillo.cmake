if(EXISTS ${ARMADILLO_INSTALL_DIR})
  message(STATUS "Armadillo is declared to be installed in ${ARMADILLO_INSTALL_DIR}")
else()
  if(ENABLE_DOWNLOAD)
    set(LIBRARIES_TO_DOWNLOAD ${LIBRARIES_TO_DOWNLOAD} ARMADILLO)
    set(ARMADILLO_URL http://sourceforge.net/projects/arma/files/armadillo-6.400.2.tar.gz)
    set(ARMADILLO_ARGS -D CMAKE_INSTALL_PREFIX=<INSTALL_DIR>)
  else(ENABLE_DOWNLOAD)
    message(FATAL_ERROR "Armadillo not found. Please provide a value to ARMADILLO_DIR or set ENABLE_DOWNLOAD to True")
  endif()
endif()
