#ifndef SIMOL_SPARSEMATRIX_HPP
#define SIMOL_SPARSEMATRIX_HPP

#include "MatrixMarketFile.hpp"

#include "eigen.hpp"

#include "dsaupd.hpp"
#include "dseupd.hpp"

namespace simol
{
  template<class ScalarType, template<class> class WrappedLibrary = eigen>
  class SparseMatrix;

  template<class ScalarType, template<class> class Library>
  std::ifstream & operator>>(std::ifstream & fileToRead, 
                             SparseMatrix<ScalarType, Library> & matrixToWrite);

  template<class ScalarType>
  class SparseMatrix<ScalarType,eigen>
  {
    friend std::ifstream & operator>> <>(std::ifstream & fileToRead, 
                                         SparseMatrix<ScalarType,eigen> & matrixToWrite);

    public:
      SparseMatrix(size_t const numberOfRows, size_t const numberOfColumns);
      SparseMatrix(MatrixMarketFile const & file);
    public:
      void eigenvalues(int numberOfEigenvaluesToCompute, 
                       std::string typeOfEigenvalues,
                       double tolerance) const;
    private:
      typename eigen<ScalarType>::SparseMatrixType wrapped_;
  };


  template<class ScalarType> inline
  SparseMatrix<ScalarType,eigen>::SparseMatrix(size_t const numberOfRows,
                                               size_t const numberOfColumns)
  : wrapped_(numberOfRows,numberOfColumns)
  {}

  template<class ScalarType> inline
  SparseMatrix<ScalarType,eigen>::SparseMatrix(MatrixMarketFile const & file)
  : wrapped_(file.numberOfRows(), file.numberOfColumns())
  {
    std::vector< Eigen::Triplet<ScalarType, std::size_t> > nonzeros(file.numberOfNonzeros());
    for(size_t nonzeroIndex = 0; nonzeroIndex < file.numberOfNonzeros(); ++nonzeroIndex)
    {
        int rowIndex;
        int columnIndex;
        ScalarType nonzero;
        fscanf(file.content(), "%d %d %lg\n", &rowIndex, &columnIndex, &nonzero);
        std::cout << rowIndex << " " << columnIndex << " " << nonzero << std::endl;
        nonzeros.push_back(Eigen::Triplet<ScalarType,size_t>(rowIndex, columnIndex, nonzero));
    }
    wrapped_.setFromTriplets(nonzeros.begin(),nonzeros.end());
  }
  
  template<class ScalarType> inline void
  SparseMatrix<ScalarType,eigen>::eigenvalues(int numberOfEigenvalues, 
                                              std::string typeOfEigenvalues,
                                              double tolerance) const
  {
    int ido, ncv, ldv, lworkl, info;
    double tol, *v;
    int *iparam, *ipntr;
    double *workd, *workl;

    ncv = 5;
    lworkl = ncv*(ncv+8);
    info = 0;
    ido = 0;
    iparam = new int[11];
    ipntr = new int[11];
    iparam[0] = 1;
    iparam[2] = 300;
    iparam[6] = 1;
    ldv = 256;
    int maxn = ldv;
    int maxncv = 25;
    v = new double[ldv*maxncv];
    workl = new double[lworkl];
    workd = new double[3*maxn];
    double * residual = new double[maxn];

    int sizeOfProblem = wrapped_.size();
    char typeOfProblem = 'I';
    dsaupd(ido, typeOfProblem, sizeOfProblem, const_cast<char*>(typeOfEigenvalues.c_str()), numberOfEigenvalues, tolerance, residual, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info);
    while( (ido==1) || (ido==-1) )
    {
      typename eigen<ScalarType>::VectorMap x(&workd[ipntr[0]-1], sizeOfProblem);
      //double *x = &workd[ipntr[0]-1];
      //double *y = &workd[ipntr[1]-1];
      auto y = wrapped_ * x;
      dsaupd(ido, typeOfProblem, sizeOfProblem, const_cast<char*>(typeOfEigenvalues.c_str()), numberOfEigenvalues, tolerance, residual, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info);
    }

    bool rvec = true;
    char howmny = 'A';
    bool *select = new bool[ncv];
    double * eigenvalues = new double[numberOfEigenvalues];
    double sigma;
    int ierr;

    dseupd( rvec, howmny, select, eigenvalues, v, ldv, sigma, typeOfProblem, sizeOfProblem, const_cast<char*>(typeOfEigenvalues.c_str()), numberOfEigenvalues, tolerance, residual, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, ierr );

    std::cout << "eigenvalues: ";
    for (int i=0;i<numberOfEigenvalues;++i)
      std::cout << eigenvalues[i] << " ";
    std::cout << std::endl; 
  }
  
  template<class ScalarType, template<class> class Library>
  std::ifstream & operator>>(std::ifstream & fileToRead, SparseMatrix<ScalarType,Library> & matrixToWrite)
  {
    std::vector< Eigen::Triplet<ScalarType,size_t> > nonzeros;

    size_t rowIndex, columnIndex;
    ScalarType value;

    while(fileToRead >> rowIndex >> columnIndex >> value)
      nonzeros.push_back(Eigen::Triplet<ScalarType,size_t>(rowIndex,columnIndex,value));

    matrixToWrite.wrapped_.setFromTriplets(nonzeros.begin(),nonzeros.end());

    return fileToRead;
  }


}


#endif
