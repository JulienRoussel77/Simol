if(EXISTS ${ARPACK_INSTALL_DIR})
  message(STATUS "Arpack is declared to be installed in ${ARPACK_INSTALL_DIR}")
else()
  if(ENABLE_DOWNLOAD)
    set(LIBRARIES_TO_DOWNLOAD ${LIBRARIES_TO_DOWNLOAD} ARPACK)
    set(ARPACK_URL http://www.caam.rice.edu/software/ARPACK/SRC/arpack96.tar.gz)
    set(ARPACK_PATCH_COMMAND cd <SOURCE_DIR> 
                          && wget http://www.caam.rice.edu/software/ARPACK/SRC/patch.tar.gz
                          && tar xvzf patch.tar.gz
                          && cp ${CMAKE_SOURCE_DIR}/cmake/arpack/CMakeLists.txt <SOURCE_DIR>)
    set(ARPACK_ARGS -D CMAKE_INSTALL_PREFIX=<INSTALL_DIR>)
    set(ARPACK_INSTALL_COMMAND mkdir -p <INSTALL_DIR>/lib && 
                               cp <BINARY_DIR>/libarpack.a <INSTALL_DIR>/lib &&
                               cp <BINARY_DIR>/libarpack_util.a <INSTALL_DIR>/lib)
  else(ENABLE_DOWNLOAD)
    message(FATAL_ERROR "Arpack not found. Please provide a value to ARPACK_INSTALL_DIR or set ENABLE_DOWNLOAD to True")
  endif()
endif()

