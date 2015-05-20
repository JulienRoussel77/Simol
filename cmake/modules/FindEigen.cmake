# Define the following variables:
# - EIGEN_INSTALL_DIR (if not user-defined)

if(EXISTS ${EIGEN_INSTALL_DIR})
  message(STATUS "Eigen is declared to be installed in ${EIGEN_INSTALL_DIR}")
else()
  if(ENABLE_DOWNLOAD)
    list(APPEND LIBRARIES_TO_DOWNLOAD EIGEN)
#    setDirectoryTree(Eigen)
    set(EIGEN_URL http://bitbucket.org/eigen/eigen/get/3.2.4.tar.gz)
    set(EIGEN_URL_MD5 4d0d77e06fef87b4fcd2c9b72cc8dc55)
    set(EIGEN_ARGS -D CMAKE_INSTALL_PREFIX=<INSTALL_DIR>)
  else(ENABLE_DOWNLOAD)
    message(FATAL_ERROR "Eigen not found. Please provide a value to EIGEN_DIR or set ENABLE_DOWNLOAD to True")
  endif()
endif()

