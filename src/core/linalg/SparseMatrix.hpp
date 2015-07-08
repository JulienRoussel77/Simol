#ifndef SIMOL_SPARSEMATRIX_HPP
#define SIMOL_SPARSEMATRIX_HPP

#include "MatrixMarketFile.hpp"


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
