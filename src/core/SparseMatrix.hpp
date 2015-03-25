#ifndef SIMOL_MATRIX_HPP
#define SIMOL_MATRIX_HPP

#include "SparseMatrix_fwd.hpp"

namespace simol
{
  template<class ScalarType, template<class> class Library>
  std::ifstream & operator>>(std::ifstream & fileToRead, SparseMatrix<ScalarType,Library> & matrixToWrite);

  template<class ScalarType>
  class SparseMatrix<ScalarType,eigen>
  {
    friend
    std::ifstream & operator>> <>(std::ifstream & fileToRead, SparseMatrix<ScalarType,eigen> & matrixToWrite);
  
    public:
      SparseMatrix(size_t const numberOfRows, size_t const numberOfColumns);
    private:
      typename eigen<ScalarType>::SparseMatrixType wrapped_;
  };

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

#include "SparseMatrix_impl.hpp"

#endif
