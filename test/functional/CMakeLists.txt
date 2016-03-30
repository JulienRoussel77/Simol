ADD_EXECUTABLE(simol_test_functional 	${CMAKE_SOURCE_DIR}/test/functional/quantchem/HartreeFockTest.cpp) 
ADD_EXECUTABLE(test_galerkin 		${CMAKE_SOURCE_DIR}/test/functional/statphys/galerkin/main.cpp) 
ADD_EXECUTABLE(test_hamiltonian 	${CMAKE_SOURCE_DIR}/test/functional/statphys/moldyn/hamiltonian/test_hamiltonian.cpp) 
ADD_EXECUTABLE(test_langevin 		${CMAKE_SOURCE_DIR}/test/functional/statphys/moldyn/langevin/test_langevin.cpp) 
ADD_EXECUTABLE(test_dpde 		${CMAKE_SOURCE_DIR}/test/functional/statphys/moldyn/dpde/test_dpde.cpp) 


INCLUDE_DIRECTORIES(	${CMAKE_SOURCE_DIR}/src/statphys/
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
TARGET_LINK_LIBRARIES(test_hamiltonian simol-statphys)
TARGET_LINK_LIBRARIES(test_langevin simol-statphys)
TARGET_LINK_LIBRARIES(test_dpde simol-statphys)

TARGET_LINK_LIBRARIES(simol_test_functional simol-quantchem)
TARGET_LINK_LIBRARIES(simol_test_functional ${EXTERNAL_DIR}/gtest/install/lib/libgtest_main.a)
TARGET_LINK_LIBRARIES(simol_test_functional ${EXTERNAL_DIR}/gtest/install/lib/libgtest.a)
TARGET_LINK_LIBRARIES(simol_test_functional ${CMAKE_THREAD_LIBS_INIT})

ENABLE_TESTING()

ADD_TEST(test_simol_functional simol_test_functional " ")

ADD_CUSTOM_TARGET(functional COMMAND simol_test_functional)
ADD_CUSTOM_TARGET(galerkin COMMAND test_galerkin -i 		"${CMAKE_SOURCE_DIR}/test/functional/statphys/galerkin/galerkin_test.yaml")
ADD_CUSTOM_TARGET(hamiltonian COMMAND test_hamiltonian -i 	"${CMAKE_SOURCE_DIR}/test/functional/statphys/moldyn/hamiltonian/hamiltonian_test.yaml")
ADD_CUSTOM_TARGET(langevin COMMAND test_langevin -i 		"${CMAKE_SOURCE_DIR}/test/functional/statphys/moldyn/langevin/langevin_sinusoidal_test.yaml")
ADD_CUSTOM_TARGET(dpde COMMAND test_dpde -i 			"${CMAKE_SOURCE_DIR}/test/functional/statphys/moldyn/dpde/dpde_1D_test.yaml")
#ADD_CUSTOM_TARGET(moldyn COMMAND test_moldyn -i tional/statphys/moldyn/langevin_sinusoidal.yaml)
