function(DownloadThenInstall LIBRARY_NAMES)

  foreach(LIBRARY_NAME IN LISTS LIBRARY_NAMES)

    string(TOLOWER ${LIBRARY_NAME} LIBRARY_LOWERNAME)
    string(TOUPPER ${LIBRARY_NAME} LIBRARY_UPPERNAME)

    ExternalProject_Add(${LIBRARY_LOWERNAME}
                        PREFIX ${EXTERNAL_DIR}/${LIBRARY_LOWERNAME}
                        DOWNLOAD_DIR ${${LIBRARY_UPPERNAME}_PREFIX}/download
                        SOURCE_DIR ${${LIBRARY_UPPERNAME}_PREFIX}/src
                        BINARY_DIR ${${LIBRARY_UPPERNAME}_PREFIX}/build
                        TMP_DIR ${${LIBRARY_UPPERNAME}_PREFIX}/tmp
                        STAMP_DIR ${${LIBRARY_UPPERNAME}_PREFIX}/stamp
                        INSTALL_DIR ${${LIBRARY_UPPERNAME}_PREFIX}/install
                        URL ${${LIBRARY_UPPERNAME}_URL}
                        URL_MD5 ${${LIBRARY_UPPERNAME}_URL_MD5}
                        CMAKE_ARGS ${${LIBRARY_UPPERNAME}_ARGS}
                        LOG_DOWNLOAD 1
                        LOG_UPDATE 1
                        LOG_CONFIGURE 1
                        LOG_BUILD 1
                        LOG_TEST 1
                        LOG_INSTALL 1
                       )
  endforeach(LIBRARY_NAME)


endfunction(DownloadThenInstall)
