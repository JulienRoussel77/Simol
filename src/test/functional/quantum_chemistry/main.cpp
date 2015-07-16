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

int main(int argc, const char* argv[])
{

  std::string const dataPath = "/home/cdoucet/Documents/soutien/matherials/codes/simol/src/test/data/quantum_chemistry/";

  simol::SparseMatrix<double> kineticMatrix(simol::MatrixMarketFile(dataPath+"kinetic_matrix.mtx"));
  simol::SparseMatrix<double> overlapMatrix(simol::MatrixMarketFile(dataPath+"overlap_matrix.mtx"));
  simol::SparseMatrix<double> potentialMatrix(simol::MatrixMarketFile(dataPath+"potential_matrix.mtx"));

  //===============
  // MATRIX LOADING
  //===============
  
/*  YAML::Node parameters = YAML::LoadFile(argv[1]);
  YAML::Node data = parameters["Data"];

  std::string path = data["path"].as<std::string>();
  
  std::ifstream file(path+data["kinetic"].as<std::string>());
  simol::SparseMatrix<double> kineticMatrix(10,10);
  file >> kineticMatrix;
  file.close();

  file.open(path+data["potential"].as<std::string>());
  simol::SparseMatrix<double> potentialMatrix(10,10);
  file >> potentialMatrix;
  file.close();

  file.open(path+data["overlap"].as<std::string>());
  simol::SparseMatrix<double> overlapMatrix(10,10);
  file >> overlapMatrix;
  file.close();

  //===========================
  // PARALLEL RANDOM GENERATION
  //===========================

  simol::Vector<double,simol::stl> vec;

  /*std::mt19937_64 mersenne;
  std::uniform_real_distribution<double> uniform(0,1);
  std::exponential_distribution<double> exponential(1);
  std::uniform_real_distribution<double> gaussian(0,1);
  #pragma omp parallel
  {
    size_t numberOfDraws = 10;
    simol::Vector<double,stl> uniform_sample(numberOfDraws);
    simol::Vector<double,stl> exponential_sample(numberOfDraws);
    simol::Vector<double,stl> gaussian_sample(numberOfDraws);
    for (size_t index = 0; index < numberOfDraws; ++index)
    {
      uniform_sample(index)     = uniform(mersenne);
      exponential_sample(index) = exponential(mersenne);
      gaussian_sample(index)    = gaussian(mersenne);
    }
    size_t identifier = omp_get_thread_num();
    #pragma omp critical
    {
      std::cout << "uniform sample of thread " << identifier << ": " << uniform_sample << std::endl;
      std::cout << "exponential sample of thread " << identifier << ": " << exponential_sample << std::endl;
      std::cout << "gaussian sample of thread " << identifier << ": " << gaussian_sample << std::endl;
      std::cout << std::endl;
    }
  }

  //==============
  // BLAS WRAPPING
  //==============

  simol::Vector<double,eigen> vec(3,1);
  vec *= 4.323;
  std::cout << vec << std::endl;

  //=======================
  // MATRIX DIAGONALIZATION
  //=======================

  simol::DenseMatrix<double> A = simol::DenseMatrix<double>::Random(6,6);
  std::cout << "Here is a random 6x6 matrix, A:" << std::endl << A << std::endl << std::endl;
  
  Eigen::EigenSolver<simol::DenseMatrix<double>> es(A);
  std::cout << "The eigenvalues of A are:" << std::endl << es.eigenvalues() << std::endl;
  std::cout << "The matrix of eigenvectors, V, is:" << std::endl << es.eigenvectors() << std::endl << std::endl;
  
  std::complex<double> lambda = es.eigenvalues()[0];
  std::cout << "Consider the first eigenvalue, lambda = " << lambda << std::endl;
  
  Eigen::VectorXcd v = es.eigenvectors().col(0);
  std::cout << "If v is the corresponding eigenvector, then lambda * v = " << std::endl << lambda * v << std::endl;
  std::cout << "... and A * v = " << std::endl << A.cast<std::complex<double> >() * v << std::endl << std::endl;
  
  Eigen::MatrixXcd D = es.eigenvalues().asDiagonal();
  Eigen::MatrixXcd V = es.eigenvectors();
  std::cout << "Finally, V * D * V^(-1) = " << std::endl << V * D * V.inverse() << std::endl;
*/
  return EXIT_SUCCESS;
}
