#include <complex>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

#include "vector.hpp"
#include "matrix.hpp"

#include "SparseMatrix.hpp"

#include <Eigen/Eigenvalues>

int main()
{

  //===============
  // MATRIX LOADING
  //===============
  
  std::string path = "/home/cdoucet/Documents/soutien/matherials/codes/existant/chimie_quantique/greedy_quantum_C++/Donnees/Donnees_H2/donnees_in";
  std::ifstream file(path+"/kinetic_matrix.txt");

  simol::SparseMatrix<double> K(10,10);
  file >> K;
  file.close();

  //===========================
  // PARALLEL RANDOM GENERATION
  //===========================

  std::mt19937_64 mersenne;
  std::uniform_real_distribution<double> uniform(0,1);
  std::exponential_distribution<double> exponential(1);
  std::uniform_real_distribution<double> gaussian(0,1);
  #pragma omp parallel
  {
    size_t numberOfDraws = 10;
    simol::Vector<double> uniform_sample(numberOfDraws);
    simol::Vector<double> exponential_sample(numberOfDraws);
    simol::Vector<double> gaussian_sample(numberOfDraws);
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

  simol::Vector<double> vec(3,1);
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

  return EXIT_SUCCESS;
}
