#ifndef SIMOL_SPARSEMATRIX_HPP
#define SIMOL_SPARSEMATRIX_HPP

namespace simol
{
  template<class ScalarType, template<class> class WrappedLibrary = eigen>
  class SparseMatrix;
}


namespace simol
{
  template<class ScalarType, template<class> class Library>
  std::ifstream & operator>>(std::ifstream & fileToRead, SparseMatrix<ScalarType,Library> & matrixToWrite);
}

namespace simol
{
  template<class ScalarType>
  class SparseMatrix<ScalarType,eigen>
  {
    friend std::ifstream & operator>> <>(std::ifstream & fileToRead, SparseMatrix<ScalarType,eigen> & matrixToWrite);

    public:
      SparseMatrix(size_t const numberOfRows, size_t const numberOfColumns);
    
    private:
      typename eigen<ScalarType>::SparseMatrixType wrapped_;
  };
}

#include "SparseMatrix.ipp"

#endif
