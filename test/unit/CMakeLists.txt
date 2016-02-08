
SIMOL_INSTALL_PACKAGE(GTest)

FIND_PACKAGE(Threads)
ADD_EXECUTABLE(simol_test_unit VectorTest.cpp)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/core/linalg)
INCLUDE_DIRECTORIES(${EXTERNAL_DIR}/gtest/install/include)

TARGET_LINK_LIBRARIES(simol_test_unit simol-core)
TARGET_LINK_LIBRARIES(simol_test_unit ${EXTERNAL_DIR}/gtest/install/lib/libgtest_main.a)
TARGET_LINK_LIBRARIES(simol_test_unit ${EXTERNAL_DIR}/gtest/install/lib/libgtest.a)
TARGET_LINK_LIBRARIES(simol_test_unit ${CMAKE_THREAD_LIBS_INIT})

