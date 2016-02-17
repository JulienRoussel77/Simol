
ADD_EXECUTABLE(simol_test_unit ${CMAKE_SOURCE_DIR}/test/unit/core/linalg/VectorTest.cpp 
                               ${CMAKE_SOURCE_DIR}/test/unit/core/linalg/DenseMatrixTest.cpp)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/core/linalg)
INCLUDE_DIRECTORIES(${EXTERNAL_DIR}/gtest/install/include)

TARGET_LINK_LIBRARIES(simol_test_unit simol-core)
TARGET_LINK_LIBRARIES(simol_test_unit ${EXTERNAL_DIR}/gtest/install/lib/libgtest_main.a)
TARGET_LINK_LIBRARIES(simol_test_unit ${EXTERNAL_DIR}/gtest/install/lib/libgtest.a)
TARGET_LINK_LIBRARIES(simol_test_unit ${CMAKE_THREAD_LIBS_INIT})

ENABLE_TESTING()

ADD_TEST(test_simol simol_test_unit " ")

ADD_CUSTOM_TARGET(unit COMMAND simol_test_unit)
