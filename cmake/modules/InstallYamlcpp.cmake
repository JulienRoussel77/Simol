if(EXISTS ${YAMLCPP_INSTALL_DIR})
  message(STATUS "Yaml-cpp is declared to be installed in ${YAMLCPP_INSTALL_DIR}")
else()
  if(ENABLE_DOWNLOAD)
    set(LIBRARIES_TO_DOWNLOAD ${LIBRARIES_TO_DOWNLOAD} YAMLCPP)
    set(YAMLCPP_URL https://yaml-cpp.googlecode.com/files/yaml-cpp-0.5.1.tar.gz)
    set(YAMLCPP_ARGS -D Boost_NO_SYSTEM_PATHS=True -D BOOST_ROOT=<INSTALL_DIR>/../../boost/install -D CMAKE_INSTALL_PREFIX=<INSTALL_DIR>)
  else(ENABLE_DOWNLOAD)
    message(FATAL_ERROR "Yaml-cpp is not found. Please provide a value to YAMLCPP_DIR or set ENABLE_DOWNLOAD to True")
  endif()
endif()




