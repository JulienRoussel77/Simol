#ifndef SIMOL_SYMSPARSEMATRIX_EIGEN_HPP
#define SIMOL_SYMSPARSEMATRIX_EIGEN_HPP

#include "simol/core/linalg/SparseMatrix.hpp"

namespace simol
{
  template<class ScalarType>
  class SymSparseMatrix<ScalarType, eigen> : public SparseMatrix<ScalarType, eigen>
  {
    public:
        SymSparseMatrix(size_t const numberOfRows, size_t const numberOfColumns);
	    virtual DenseMatrix<ScalarType, eigen> dense() const;
  };

    //! Returns a dense matrix equal to the sparse matrix
    template<typename ScalarType> inline
    DenseMatrix<ScalarType, eigen> SymSparseMatrix<ScalarType, eigen>::dense() const
    {
       DenseMatrix<ScalarType, eigen> M(this->wrapped_.rows(), this->wrapped_.cols());
       M.wrapped_ = this->wrapped_.template triangularView<Eigen::Upper>();
       M.wrapped_ += (this->wrapped_.transpose()).template triangularView<Eigen::StrictlyLower>();
       return M;
    }

    //! Construction from a given size
    template<class ScalarType> inline
    SymSparseMatrix<ScalarType, eigen>::SymSparseMatrix(size_t const numberOfRows, size_t const numberOfColumns)
      : SparseMatrix<ScalarType, eigen>(numberOfRows, numberOfColumns)
    {}


}


#endif

