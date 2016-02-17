
ADD_EXECUTABLE(simol_test_quantchem ${CMAKE_SOURCE_DIR}/test/functional/quantchem/HartreeFockTest.cpp)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/quantchem)
INCLUDE_DIRECTORIES(${EXTERNAL_DIR}/gtest/install/include)

TARGET_LINK_LIBRARIES(simol_test_quantchem simol-quantchem)
TARGET_LINK_LIBRARIES(simol_test_quantchem ${EXTERNAL_DIR}/gtest/install/lib/libgtest_main.a)
TARGET_LINK_LIBRARIES(simol_test_quantchem ${EXTERNAL_DIR}/gtest/install/lib/libgtest.a)
TARGET_LINK_LIBRARIES(simol_test_quantchem ${CMAKE_THREAD_LIBS_INIT})

ENABLE_TESTING()

ADD_TEST(test_simol_quantchem simol_test_quantchem " ")

add_custom_target(test_quantchem COMMAND simol_test_quantchem)

