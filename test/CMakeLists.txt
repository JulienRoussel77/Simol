SIMOL_INSTALL_PACKAGE(GTest)

FIND_PACKAGE(Threads)
#add_executable(test_simol unit/test_simol.cpp 
#                          unit/EigenSolverTest.cpp 
#                          unit/SparseMatrixTest.cpp 
#                          ../src/core/io/mmio.c
#                          unit/VectorTest.cpp)

#target_link_libraries(test_simol ${CPPUNIT_INSTALL_DIR}/lib/libcppunit.so 
#                                 ${ARPACK_INSTALL_DIR}/lib/libarpack.a
#                                 ${ARPACK_INSTALL_DIR}/lib/libarpack_util.a
#                                 ${LAPACK_LIB}
#                                 ${BLAS_LIB}
#                                 )
ENABLE_TESTING()

#add_test(test_simol ${EXECUTABLE_OUTPUT_PATH}/test_simol)

#set_tests_properties(test_simol PROPERTIES PASS_REGULAR_EXPRESSION "OK")

ADD_SUBDIRECTORY(unit)
ADD_SUBDIRECTORY(functional)
