#include <complex>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

#include "linalg/SparseMatrix.hpp"
#include <Eigen/Eigenvalues>

#include <yaml-cpp/yaml.h>
#include "io/MatrixMarketFile.hpp"

#include "../../../core/quantchem/hartree_fock.hpp"

int main()
{

  std::string const dataPath = "/home/cdoucet/Documents/soutien/matherials/codes/simol/src/test/data/quantum_chemistry/";

  simol::SparseMatrix<double> kineticMatrix(simol::MatrixMarketFile(dataPath+"kinetic_matrix.mtx"));
  simol::SparseMatrix<double> overlapMatrix(simol::MatrixMarketFile(dataPath+"overlap_matrix.mtx"));
  simol::SparseMatrix<double> potentialMatrix(simol::MatrixMarketFile(dataPath+"potential_matrix.mtx"));
  return EXIT_SUCCESS;
}

