if(EXISTS ${BOOST_INSTALL_DIR})
  message(STATUS "Boost is declared to be installed in ${BOOST_INSTALL_DIR}")
else()
  if(ENABLE_DOWNLOAD)
    set(LIBRARIES_TO_DOWNLOAD ${LIBRARIES_TO_DOWNLOAD} BOOST)
    set(BOOST_URL http://sourceforge.net/projects/boost/files/boost/1.57.0/boost_1_57_0.tar.gz)
    set(BOOST_URL_MD5 25f9a8ac28beeb5ab84aa98510305299)
    set(BOOST_CONFIGURE_COMMAND pwd)
    set(BOOST_BUILD_COMMAND pwd)
    set(BOOST_INSTALL_COMMAND mkdir <INSTALL_DIR>/include && cp -r <SOURCE_DIR>/boost <INSTALL_DIR>/include/boost)
  else(ENABLE_DOWNLOAD)
    message(FATAL_ERROR "Boost not found. Please provide a value to BOOST_DIR or set ENABLE_DOWNLOAD to True")
  endif()
endif()

