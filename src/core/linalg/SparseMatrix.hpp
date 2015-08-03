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
      std::size_t numberOfRows() const;
      std::size_t numberOfColumns() const;
    public:
      typename eigen<ScalarType>::SparseMatrixType const & wrapped() const
      { return wrapped_; }
    public:
      SparseMatrix(size_t const numberOfRows, size_t const numberOfColumns);
      SparseMatrix(MatrixMarketFile const & file);
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
        nonzeros.push_back(Eigen::Triplet<ScalarType,size_t>(rowIndex, columnIndex, nonzero));
    }
    wrapped_.setFromTriplets(nonzeros.begin(),nonzeros.end());
  }

  template<class ScalarType> inline std::size_t
  SparseMatrix<ScalarType,eigen>::numberOfRows() const
  { return wrapped_.rows(); }

  template<class ScalarType> inline std::size_t
  SparseMatrix<ScalarType,eigen>::numberOfColumns() const
  { return wrapped_.cols(); }

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
