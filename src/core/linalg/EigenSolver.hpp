#ifndef SIMOL_EIGENSOLVER_HPP
#define SIMOL_EIGENSOLVER_HPP

#include "SparseMatrix.hpp"


namespace simol
{

  template<typename ScalarType>
  struct arpack
  {};

  template<typename ScalarType, template<class> class WrappedEigenSolver = arpack>
  class EigenSolver
  {
    public:
      double* solve(SparseMatrix<ScalarType, eigen> const & A, 
                    int numberOfEigenvalues,
                    double tolerance) const;
  };

  template<typename ScalarType, template<typename> class WrappedEigenSolver> double*
  EigenSolver<ScalarType, WrappedEigenSolver>::solve(SparseMatrix<ScalarType, eigen> const & leftMatrix, 
                                                     int numberOfEigenvalues,
                                                     double tolerance) const
  {

    int problemSize = leftMatrix.numberOfRows();
    int ncv = 5;
    char bmat = 'I';
    char * whichEigenvalues = "LM";
    int info = 0;
  
    int iparam[11];
    iparam[0] = 1;
    iparam[2] = 300;
    iparam[6] = 1;
    
    int ipntr[11];
    
    int ldv = 256;
    int maxn = ldv;
    int maxncv = 25;
    double * v = new double[ldv*maxncv];
    
    // The length of workl is always given by the above formula (see Arpack Users' guide)
    int lworkl = ncv*(ncv+8);
    double * workl = new double[lworkl];
    
    double * workd = new double[3*maxn];
    double * residual = new double[maxn];

    int ido = 0;
    dsaupd(ido, bmat, problemSize, whichEigenvalues, numberOfEigenvalues, tolerance, residual, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info);
    while( (ido==1) || (ido==-1) )
    {
      typename eigen<ScalarType>::VectorMap x(&workd[ipntr[0]-1], problemSize);
      typename eigen<ScalarType>::VectorMap y(&workd[ipntr[1]-1], problemSize);
      y = leftMatrix.wrapped() * x;
      dsaupd(ido, bmat, problemSize, whichEigenvalues, numberOfEigenvalues, tolerance, residual, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info);

    }

    bool rvec = true;
    char howmny = 'A';
    bool *select = new bool[ncv];
    double * eigenvalues = new double[numberOfEigenvalues];
    double sigma;
    int ierr;

    dseupd( rvec, howmny, select, eigenvalues, v, ldv, sigma, bmat, problemSize, whichEigenvalues, numberOfEigenvalues, tolerance, residual, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, ierr );
    
    return eigenvalues;

  }

}

#endif
