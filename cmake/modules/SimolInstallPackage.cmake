function(simol_install_package LIBRARY_NAMES)

foreach(LIBRARY_NAME IN LISTS LIBRARY_NAMES)

    include(Install${LIBRARY_NAME})

    string(TOLOWER ${LIBRARY_NAME} LIBRARY_LOWERNAME)
    string(TOUPPER ${LIBRARY_NAME} LIBRARY_UPPERNAME)

    set(${LIBRARY_UPPERNAME}_PREFIX ${SIMOL_EXTERNAL_DIR}/${LIBRARY_LOWERNAME} )
    set(${LIBRARY_UPPERNAME}_DOWNLOAD_DIR ${${LIBRARY_UPPERNAME}_PREFIX}/download)
    set(${LIBRARY_UPPERNAME}_SOURCE_DIR ${${LIBRARY_UPPERNAME}_PREFIX}/src)
    set(${LIBRARY_UPPERNAME}_BINARY_DIR ${${LIBRARY_UPPERNAME}_PREFIX}/build)
    set(${LIBRARY_UPPERNAME}_TMP_DIR ${${LIBRARY_UPPERNAME}_PREFIX}/tmp)
    set(${LIBRARY_UPPERNAME}_STAMP_DIR ${${LIBRARY_UPPERNAME}_PREFIX}/stamp)
    set(${LIBRARY_UPPERNAME}_INSTALL_DIR ${${LIBRARY_UPPERNAME}_PREFIX}/install)

    ExternalProject_Add(${LIBRARY_LOWERNAME}
                        PREFIX ${${LIBRARY_UPPERNAME}_PREFIX}
                        DOWNLOAD_DIR ${${LIBRARY_UPPERNAME}_DOWNLOAD_DIR}
                        SOURCE_DIR ${${LIBRARY_UPPERNAME}_SOURCE_DIR}
                        BINARY_DIR ${${LIBRARY_UPPERNAME}_BINARY_DIR}
                        TMP_DIR ${${LIBRARY_UPPERNAME}_TMP_DIR}
                        STAMP_DIR ${${LIBRARY_UPPERNAME}_STAMP_DIR}
                        INSTALL_DIR ${${LIBRARY_UPPERNAME}_INSTALL_DIR}
                        URL ${${LIBRARY_UPPERNAME}_URL}
                        CMAKE_ARGS ${${LIBRARY_UPPERNAME}_ARGS}
                        PATCH_COMMAND ${${LIBRARY_UPPERNAME}_PATCH_COMMAND}
                        CONFIGURE_COMMAND ${${LIBRARY_UPPERNAME}_CONFIGURE_COMMAND}
                        BUILD_COMMAND ${${LIBRARY_UPPERNAME}_BUILD_COMMAND}
                        INSTALL_COMMAND ${${LIBRARY_UPPERNAME}_INSTALL_COMMAND}
               #         LOG_DOWNLOAD 1
                        LOG_UPDATE 1
               #         LOG_CONFIGURE 1
               #         LOG_BUILD 1
                        LOG_TEST 1
              #          LOG_INSTALL 1
                       )
    set(${LIBRARY_UPPERNAME}_INSTALL_DIR ${${LIBRARY_UPPERNAME}_PREFIX}/install PARENT_SCOPE)
    message(STATUS "${LIBRARY_UPPERNAME}_INSTALL_DIR is ${${LIBRARY_UPPERNAME}_INSTALL_DIR}")

endforeach(LIBRARY_NAME)


endfunction(simol_install_package)
