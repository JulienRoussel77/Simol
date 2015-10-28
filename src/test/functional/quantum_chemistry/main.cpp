#include <complex>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

#include "linalg/SparseMatrix.hpp"
#include "quantchem/SparseTensor.hpp"
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
  
  std::string filename = dataPath + "n.txt";
  FILE* size_file = fopen(filename.c_str(), "r");
  size_t M_disc;
  fscanf(size_file, "%lu", &M_disc);
  
  simol::DenseMatrix<double> kinetic_matrix(kinetic_file);
  simol::DenseMatrix<double> overlap_matrix(overlap_file);
  simol::DenseMatrix<double> potential_matrix(potential_file);
  simol::SparseTensor<double> twoelectron_matrix(dataPath + "twoelectron_matrix.txt", M_disc);
  
  /*simol::SparseMatrix<double> kineticMatrix(simol::MatrixMarketFile(dataPath+"kinetic_matrix.mtx"));
  simol::SparseMatrix<double> overlapMatrix(simol::MatrixMarketFile(dataPath+"overlap_matrix.mtx"));
  simol::SparseMatrix<double> potentialMatrix(simol::MatrixMarketFile(dataPath+"potential_matrix.mtx"));*/
  return EXIT_SUCCESS;
}

