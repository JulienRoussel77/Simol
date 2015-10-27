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

  std::string const dataPath = "/home/cdoucet/Documents/soutien/matherials/simol/src/test/data/LiH/";

  simol::MatrixMarketFile kinetic_file(dataPath + "kinetic_matrix.mtx");
  simol::MatrixMarketFile overlap_file(dataPath + "overlap_matrix.mtx");
  simol::MatrixMarketFile potential_file(dataPath + "potential_matrix.mtx");
  simol::MatrixMarketFile twoelectron_file(dataPath + "twoelectron_matrix.mtx");
  
  simol::DenseMatrix<double> kinetic_matrix(kinetic_file);
  simol::DenseMatrix<double> overlap_matrix(overlap_file);
  simol::DenseMatrix<double> potential_matrix(potential_file);
  simol::SparseMatrix<double> twoelectron_matrix(twoelectron_file);
  
  /*simol::SparseMatrix<double> kineticMatrix(simol::MatrixMarketFile(dataPath+"kinetic_matrix.mtx"));
  simol::SparseMatrix<double> overlapMatrix(simol::MatrixMarketFile(dataPath+"overlap_matrix.mtx"));
  simol::SparseMatrix<double> potentialMatrix(simol::MatrixMarketFile(dataPath+"potential_matrix.mtx"));*/
  return EXIT_SUCCESS;
}

