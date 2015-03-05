#include <complex>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

#include "core/Vector.hpp"
#include "core/matrix.hpp"

#include "core/SparseMatrix.hpp"
#include <Eigen/Eigenvalues>

#include <yaml-cpp/yaml.h>

extern "C"
{

  void dsaupd(int *ido, char *bmat, int *n, char *which,
                       int *nev, double *tol, double *resid,
                       int *ncv, double *v, int *ldv,
                       int *iparam, int *ipntr, double *workd,
                       double *workl, int *lworkl, int *info);
}




int main(int argc, const char* argv[])
{


/*
  //=======
  // ARPACK
  //=======
  int ido = 0;
  char bmat = 'I';
  int n = 3;
  char* which = "LM";
  int nev = 1;
  double tol = 1e-6;
  double* resid;
  int ncv;
  double* v;
  int ldv;
  int* iparam;
  int* ipntr;
  double* workd;
  double* workl;
  int lworkl;
  int info;
*/
  /*dsapud(ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info);

 // a.1) User defined parameters.

  int     n;          // Dimension of the eigenproblem.
  int     nev;        // Number of eigenvalues to be computed. 0 < nev < n-1.
  int     ncv;        // Number of Arnoldi vectors generated at each iteration.
  int     maxit;      // Maximum number of Arnoldi update iterations allowed.
  char*   which;      // Specify which of the Ritz values of OP to compute.
  FLOAT   tol;        // Stopping criterion (relative accuracy of Ritz values).
  FLOAT   sigmaI;     // Imaginary part of shift (for nonsymmetric problems).
  TYPE    sigmaR;     // Shift (real part only if problem is nonsymmetric).
  TYPE    *resid;     // Initial residual vector.


 // a.2) Internal variables.

  bool    rvec;       // Indicates if eigenvectors/Schur vectors were
                      // requested (or only eigenvalues will be determined).
  bool    newRes;     // Indicates if a new "resid" vector was created.
  bool    newVal;     // Indicates if a new "EigValR" vector was created.
  bool    newVec;     // Indicates if a new "EigVec" vector was created.
  bool    PrepareOK;  // Indicates if internal variables were correctly set.
  bool    BasisOK;    // Indicates if an Arnoldi basis was found.
  bool    ValuesOK;   // Indicates if eigenvalues were calculated.
  bool    VectorsOK;  // Indicates if eigenvectors were determined.
  bool    SchurOK;    // Indicates if Schur vectors were determined.
  bool    AutoShift;  // Indicates if implicit shifts will be generated
                      // internally (or will be supplied by the user).
  char    bmat;       // Indicates if the problem is a standard ('I') or
                      // generalized ('G") eigenproblem.
  char    HowMny;     // Indicates if eigenvectors ('A') or Schur vectors ('P')
                      // were requested (not referenced if rvec = false).
  int     ido;        // Original ARPACK reverse communication flag.
  int     info;       // Original ARPACK error flag.
  int     mode;       // Indicates the type of the eigenproblem (regular,
                      // shift and invert, etc).
  int     lworkl;     // Dimension of array workl.
  int     lworkv;     // Dimension of array workv.
  int     lrwork;     // Dimension of array rwork.
  int     iparam[12]; // Vector that handles original ARPACK parameters.
  int     ipntr[15];  // Vector that handles original ARPACK pointers.
  FLOAT   *rwork;     // Original ARPACK internal vector.
  TYPE    *workl;     // Original ARPACK internal vector.
  TYPE    *workd;     // Original ARPACK internal vector.
  TYPE    *workv;     // Original ARPACK internal vector.
  TYPE    *V;         // Arnoldi basis / Schur vectors.*/

  //===============
  // MATRIX LOADING
  //===============
  
  std::ifstream file(argv[1]);
  simol::SparseMatrix<double,simol::fromEigen> kineticMatrix(10,10);
  file >> kineticMatrix;
  file.close();

  file.open(argv[2]);
  simol::SparseMatrix<double,simol::fromEigen> potentialMatrix(10,10);
  file >> potentialMatrix;
  file.close();

  file.open(argv[3]);
  simol::SparseMatrix<double,simol::fromEigen> overlapMatrix(10,10);
  file >> overlapMatrix;
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
    simol::Vector<double,simol::Dummy> uniform_sample(numberOfDraws);
    simol::Vector<double,simol::Dummy> exponential_sample(numberOfDraws);
    simol::Vector<double,simol::Dummy> gaussian_sample(numberOfDraws);
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

  simol::Vector<double,simol::Dummy> vec(3,1);
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
