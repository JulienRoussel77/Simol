
if(EXISTS ${GTEST_INSTALL_DIR})
  message(STATUS "Google test is declared to be installed in ${GTEST_INSTALL_DIR}")
else()
  if(ENABLE_DOWNLOAD)
    set(LIBRARIES_TO_DOWNLOAD ${LIBRARIES_TO_DOWNLOAD} GTEST)
    set(GTEST_URL https://github.com/google/googletest/archive/release-1.7.0.tar.gz)
    set(GTEST_INSTALL_COMMAND mkdir -p ../install/lib && cp libgtest.a ../install/lib && cp libgtest_main.a ../install/lib && cp -r ../src/include ../install )
  else(ENABLE_DOWNLOAD)
      message(FATAL_ERROR "Google Test (gtest) not found. Please provide a value to GTEST_DIR or set ENABLE_DOWNLOAD to True")
  endif()
endif()

