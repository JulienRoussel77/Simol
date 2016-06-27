#ifndef SIMOL_SYMSPARSEMATRIX_EIGEN_HPP
#define SIMOL_SYMSPARSEMATRIX_EIGEN_HPP

#include "simol/core/linalg/SparseMatrix.hpp"

namespace simol
{
  template<class ScalarType>
  class SymSparseMatrix<ScalarType, eigen> : public SparseMatrix<ScalarType, eigen>
  {
    public:
      explicit SymSparseMatrix(size_t const numberOfRows, size_t const numberOfColumns);
      explicit SymSparseMatrix(MatrixMarketFile const & file);
      explicit SymSparseMatrix(std::string const & filename, std::size_t const size);
      //explicit SymSparseMatrix(Vector<double, eigen> u, std::size_t numberOfRows, std::size_t numberOfColumns);
	
        
        virtual DenseMatrix<ScalarType, eigen> dense() const;
  };

    //! Construction from a given size
    template<class ScalarType> inline
    SymSparseMatrix<ScalarType, eigen>::SymSparseMatrix(size_t const numberOfRows, size_t const numberOfColumns)
      : SparseMatrix<ScalarType, eigen>(numberOfRows, numberOfColumns)
    {}
    
    //! Construction from a matrix Market file
    template<class ScalarType> inline
    SymSparseMatrix<ScalarType, eigen>::SymSparseMatrix(MatrixMarketFile const & file)
      : SparseMatrix<ScalarType, eigen>(file)
    {}
    
    //! Construction from a file
    template<class ScalarType> inline
    SymSparseMatrix<ScalarType, eigen>::SymSparseMatrix(std::string const & filename, std::size_t const size)
      : SparseMatrix<ScalarType, eigen>(filename,size)
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

