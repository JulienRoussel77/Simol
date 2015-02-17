#ifndef SPARSEMATRIX_HPP
#define SPRASEMATRIX_HPP

#include<Eigen/Sparse>
#include<Eigen/SparseCore>

#include<fstream>

#include "SparseMatrix_fwd.hpp"

namespace simol
{

  template<class ScalarType>
  std::ifstream & operator>>(std::ifstream & fileToRead, SparseMatrix<ScalarType> & matrixToWrite);

  template<class ScalarType>
  class SparseMatrix
  {
    friend std::ifstream & operator>> <>(std::ifstream &, SparseMatrix<ScalarType> &);

  public:

    SparseMatrix() = delete;
    SparseMatrix(size_t numberOfRows, size_t numberOfColumns);
    SparseMatrix(SparseMatrix const &) = delete;
    SparseMatrix & operator=(SparseMatrix const &) = delete;
    ~SparseMatrix() = default;

  private:

    Eigen::SparseMatrix<ScalarType> wrapped_;
  };

  #include "SparseMatrix_inline.hpp"

  template<class ScalarType>
  std::ifstream & operator>>(std::ifstream & fileToRead, SparseMatrix<ScalarType> & matrixToWrite)
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
