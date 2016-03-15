ADD_EXECUTABLE(simol_test_functional ${CMAKE_SOURCE_DIR}/test/functional/quantchem/HartreeFockTest.cpp) 
ADD_EXECUTABLE(test_galerkin ${CMAKE_SOURCE_DIR}/test/functional/statphys/galerkin/main.cpp) 
ADD_EXECUTABLE(test_moldyn ${CMAKE_SOURCE_DIR}/test/functional/statphys/moldyn/main.cpp) 



INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/statphys/
                    ${CMAKE_SOURCE_DIR}/src/statphys/input
		${CMAKE_SOURCE_DIR}/src/statphys/output
		${CMAKE_SOURCE_DIR}/src/statphys/dynamics
		${CMAKE_SOURCE_DIR}/src/statphys/system	
		${CMAKE_SOURCE_DIR}/src/statphys/controlVariate
		${CMAKE_SOURCE_DIR}/src/statphys/potential
		${CMAKE_SOURCE_DIR}/src/core/linalg
		${CMAKE_SOURCE_DIR}/src/core/random
		${CMAKE_SOURCE_DIR}/src/core/io
		${CMAKE_SOURCE_DIR}/src/core
)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/quantchem)
INCLUDE_DIRECTORIES(${EXTERNAL_DIR}/gtest/install/include)


TARGET_LINK_LIBRARIES(test_galerkin simol-statphys)
TARGET_LINK_LIBRARIES(test_moldyn simol-statphys)

TARGET_LINK_LIBRARIES(simol_test_functional simol-quantchem)
TARGET_LINK_LIBRARIES(simol_test_functional ${EXTERNAL_DIR}/gtest/install/lib/libgtest_main.a)
TARGET_LINK_LIBRARIES(simol_test_functional ${EXTERNAL_DIR}/gtest/install/lib/libgtest.a)
TARGET_LINK_LIBRARIES(simol_test_functional ${CMAKE_THREAD_LIBS_INIT})

ENABLE_TESTING()

ADD_TEST(test_simol_functional simol_test_functional " ")

ADD_CUSTOM_TARGET(functional COMMAND simol_test_functional)
ADD_CUSTOM_TARGET(galerkin COMMAND test_galerkin -i "${CMAKE_SOURCE_DIR}/test/functional/statphys/galerkin/galerkin_test.yaml")
ADD_CUSTOM_TARGET(moldyn COMMAND test_moldyn -i "${CMAKE_SOURCE_DIR}/test/functional/statphys/moldyn/langevin_sinusoidal_test.yaml")
#ADD_CUSTOM_TARGET(moldyn COMMAND test_moldyn -i tional/statphys/moldyn/langevin_sinusoidal.yaml)
