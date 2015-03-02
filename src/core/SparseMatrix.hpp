#ifndef SPARSEMATRIX_HPP
#define SPRASEMATRIX_HPP


#include<fstream>

#include "../fromEigen.hpp"
#include "SparseMatrix_fwd.hpp"

namespace simol
{
  template<class ScalarType, template<class> class WrappingPolicy>
  std::ifstream & operator>>(std::ifstream & fileToRead, SparseMatrix<ScalarType,WrappingPolicy> & matrixToWrite);
}

namespace simol
{

  template<class ScalarType, template<class> class WrappingPolicy>
  class SparseMatrix : public WrappingPolicy<ScalarType>
  {
    
    //=================
    // FRIEND FUNCTIONS
    //=================

    friend std::ifstream & operator>> <>(std::ifstream & fileToRead, SparseMatrix<ScalarType,WrappingPolicy> & matrixToWrite);

  public:

    SparseMatrix() = delete;
    SparseMatrix(size_t numberOfRows, size_t numberOfColumns);
    SparseMatrix(SparseMatrix const &) = delete;
    SparseMatrix & operator=(SparseMatrix const &) = delete;
    ~SparseMatrix() = default;

  private:

    typename WrappingPolicy<ScalarType>::WrappedType wrapped_;
  };

  #include "SparseMatrix_impl.hpp"

template<class ScalarType, template<class> class WrappingPolicy>
std::ifstream & operator>>(std::ifstream & fileToRead, SparseMatrix<ScalarType,WrappingPolicy> & matrixToWrite)
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
