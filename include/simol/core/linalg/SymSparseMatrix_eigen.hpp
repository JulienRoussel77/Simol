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

    //! Construction from a given size
    template<class ScalarType> inline
    SymSparseMatrix<ScalarType, eigen>::SymSparseMatrix(size_t const numberOfRows, size_t const numberOfColumns)
      : SparseMatrix<ScalarType, eigen>(numberOfRows, numberOfColumns)
    {}
  
  
    //! Returns a dense matrix equal to the sparse matrix
    template<typename ScalarType> inline
    DenseMatrix<ScalarType, eigen> SymSparseMatrix<ScalarType, eigen>::dense() const
    {
       DenseMatrix<ScalarType, eigen> M(this->wrapped_.rows(), this->wrapped_.cols());
       M.wrapped_ = this->wrapped_.template triangularView<Eigen::Upper>();
       M.wrapped_ += (this->wrapped_.transpose()).template triangularView<Eigen::StrictlyLower>();
       return M;
    }
    
    //! Do a matrix by vector product: symmetric version 
    template<class ScalarType>
    Vector<ScalarType, eigen> operator*(SymSparseMatrix<ScalarType, eigen> matrix, Vector<ScalarType, eigen> const & vector)
    {
        Vector<ScalarType, eigen> prod(matrix.wrapped_.rows());
        prod.wrapped_ = matrix.wrapped_.template selfadjointView<Eigen::Upper>() * vector.wrapped_;
        return prod;
     }
}


#endif

