#ifndef SIMOL_SPARSEMATRIX_IMPL_HPP
#define SIMOL_SPARSEMATRIX_IMPL_HPP

#include <cstdlib>


namespace simol
{

  //=============
  // CONSTRUCTORS
  //=============

  template<class ScalarType> inline
  SparseMatrix<ScalarType,eigen>::SparseMatrix(size_t const numberOfRows,
                                               size_t const numberOfColumns)
  :wrapped_(numberOfRows,numberOfColumns)
  {}

  //=============
  // INPUT/OUTPUT
  //=============
  
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
