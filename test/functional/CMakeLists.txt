ADD_EXECUTABLE(simol_test_functional 	${CMAKE_SOURCE_DIR}/test/functional/quantchem/HartreeFockTest.cpp) 
ADD_EXECUTABLE(test_galerkin 		${CMAKE_SOURCE_DIR}/test/functional/statphys/galerkin/TestGalerkin.cpp) 
ADD_EXECUTABLE(test_hamiltonian 	${CMAKE_SOURCE_DIR}/test/functional/statphys/moldyn/hamiltonian/TestHamiltonian.cpp) 
ADD_EXECUTABLE(test_langevin 		${CMAKE_SOURCE_DIR}/test/functional/statphys/moldyn/langevin/TestLangevin.cpp) 
ADD_EXECUTABLE(test_dpde 		${CMAKE_SOURCE_DIR}/test/functional/statphys/moldyn/dpde/TestDPDE.cpp) 
ADD_EXECUTABLE(test_rotor 		${CMAKE_SOURCE_DIR}/test/functional/statphys/moldyn/chain/TestRotor.cpp)
ADD_EXECUTABLE(test_bichainfpu 		${CMAKE_SOURCE_DIR}/test/functional/statphys/moldyn/chain/TestBiChainFPU.cpp)
ADD_EXECUTABLE(test_trichain 		${CMAKE_SOURCE_DIR}/test/functional/statphys/moldyn/chain/TestTriChain.cpp)
ADD_EXECUTABLE(test_fluid 		${CMAKE_SOURCE_DIR}/test/functional/statphys/moldyn/fluid/TestFluid.cpp)  

INCLUDE_DIRECTORIES(	${CMAKE_SOURCE_DIR}/src/statphys/
	    		${CMAKE_SOURCE_DIR}/src/statphys/input
			${CMAKE_SOURCE_DIR}/src/statphys/output
			${CMAKE_SOURCE_DIR}/src/statphys/dynamics
			${CMAKE_SOURCE_DIR}/src/statphys/system	
			${CMAKE_SOURCE_DIR}/src/statphys/simulation	
			${CMAKE_SOURCE_DIR}/src/statphys/controlVariate
			${CMAKE_SOURCE_DIR}/src/statphys/potential
			${CMAKE_SOURCE_DIR}/src/core/linalg
			${CMAKE_SOURCE_DIR}/src/core/random
			${CMAKE_SOURCE_DIR}/src/core/io
			${CMAKE_SOURCE_DIR}/src/core
)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/quantchem)
INCLUDE_DIRECTORIES(${SIMOL_EXTERNAL_DIR}/gtest/install/include)


TARGET_LINK_LIBRARIES(test_galerkin simol-statphys)
TARGET_LINK_LIBRARIES(test_hamiltonian simol-statphys)
TARGET_LINK_LIBRARIES(test_langevin simol-statphys)
TARGET_LINK_LIBRARIES(test_dpde simol-statphys)
TARGET_LINK_LIBRARIES(test_rotor simol-statphys)
TARGET_LINK_LIBRARIES(test_bichainfpu simol-statphys)
TARGET_LINK_LIBRARIES(test_trichain simol-statphys)
TARGET_LINK_LIBRARIES(test_fluid simol-statphys)

TARGET_LINK_LIBRARIES(simol_test_functional simol-quantchem)
TARGET_LINK_LIBRARIES(simol_test_functional ${SIMOL_EXTERNAL_DIR}/gtest/install/lib/libgtest_main.a)
TARGET_LINK_LIBRARIES(simol_test_functional ${SIMOL_EXTERNAL_DIR}/gtest/install/lib/libgtest.a)
TARGET_LINK_LIBRARIES(simol_test_functional ${CMAKE_THREAD_LIBS_INIT})

ENABLE_TESTING()

ADD_TEST(test_simol_functional simol_test_functional " ")

ADD_CUSTOM_TARGET(functional COMMAND simol_test_functional)
ADD_CUSTOM_TARGET(galerkin COMMAND test_galerkin -i 		"${CMAKE_SOURCE_DIR}/test/functional/statphys/galerkin/galerkin_test.yaml")
ADD_CUSTOM_TARGET(hamiltonian COMMAND test_hamiltonian -i 	"${CMAKE_SOURCE_DIR}/test/functional/statphys/moldyn/hamiltonian/hamiltonian_test.yaml")
ADD_CUSTOM_TARGET(langevin COMMAND test_langevin -i 		"${CMAKE_SOURCE_DIR}/test/functional/statphys/moldyn/langevin/langevin_sinusoidal_test.yaml")
ADD_CUSTOM_TARGET(dpde COMMAND test_dpde -i 			"${CMAKE_SOURCE_DIR}/test/functional/statphys/moldyn/dpde/dpde_1D_test.yaml")
ADD_CUSTOM_TARGET(rotor COMMAND test_rotor -i 			"${CMAKE_SOURCE_DIR}/test/functional/statphys/moldyn/chain/rotor_test.yaml")
ADD_CUSTOM_TARGET(bichainfpu COMMAND test_bichainfpu -i 			"${CMAKE_SOURCE_DIR}/test/functional/statphys/moldyn/chain/bichainfpu_test.yaml")
ADD_CUSTOM_TARGET(trichain COMMAND test_trichain -i 		"${CMAKE_SOURCE_DIR}/test/functional/statphys/moldyn/chain/trichain_test.yaml")
ADD_CUSTOM_TARGET(fluid COMMAND test_fluid -i 			"${CMAKE_SOURCE_DIR}/test/functional/statphys/moldyn/fluid/fluid_test.yaml")
