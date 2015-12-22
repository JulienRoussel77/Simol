if(EXISTS ${CPPUNIT_INSTALL_DIR})
  message(STATUS "CppUnit is declared to be installed in ${CPPUNIT_INSTALL_DIR}")
else()
  if(ENABLE_DOWNLOAD)
    set(LIBRARIES_TO_DOWNLOAD ${LIBRARIES_TO_DOWNLOAD} CPPUNIT)
    set(CPPUNIT_URL http://sourceforge.net/projects/cppunit/files/cppunit/1.12.1/cppunit-1.12.1.tar.gz)
    set(CPPUNIT_CONFIGURE_COMMAND LDFLAGS=-ldl CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} <SOURCE_DIR>/configure --prefix=<INSTALL_DIR>)
 #   set(EIGEN_ARGS -D CMAKE_INSTALL_PREFIX=<INSTALL_DIR>)
  else(ENABLE_DOWNLOAD)
    message(FATAL_ERROR "Eigen not found. Please provide a value to EIGEN_DIR or set ENABLE_DOWNLOAD to True")
  endif()
endif()
#                      URL 
#                      CONFIGURE_COMMAND cd ../cppunit && LDFLAGS="-ldl" ./configure --prefix=${PREFIX}

