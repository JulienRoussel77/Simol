#ifndef SIMOL_DENSEMATRIX_HPP
#define SIMOL_DENSEMATRIX_HPP

#include "simol/core/io/MatrixMarketFile.hpp"
#include "simol/core/linalg/Vector.hpp"

#include <ostream>

namespace simol
{

  template<typename Scalar, typename Wrapper>
  class SparseMatrix;
  
  template<typename Scalar, typename Wrapper = eigen>
  class DenseMatrix
  {
    public:
      
      static DenseMatrix Zero(std::size_t const numberOfRows, std::size_t numberOfColumns);
      static DenseMatrix Identity(std::size_t const dimension);

      DenseMatrix(DenseMatrix const & matrix) = default;
      ~DenseMatrix() = default;
      
      DenseMatrix(MatrixMarketFile const & file);
      DenseMatrix(SparseMatrix<Scalar, Wrapper> const & matrix);
      DenseMatrix(std::size_t numberOfRows, std::size_t numberOfColumns);
      DenseMatrix(typename Wrapper::template DenseMatrix<Scalar> const & matrix);
      DenseMatrix(typename Wrapper::template DenseBlock_const<Scalar> const & block);
      DenseMatrix(Vector<Scalar, Wrapper> u, std::size_t const numberOfRows, std::size_t const numberOfColumns);

      void fill(Scalar const scalar);

      std::size_t number_of_rows() const;
      std::size_t number_of_columns() const;
      
      Scalar determinant() const;
      Scalar rcond() const;
      Scalar trace() const;

      Scalar & operator()(size_t const rowIndex, size_t const columnIndex);
      
      Scalar const & operator()(size_t const rowIndex, size_t const columnIndex) const;

      Vector<Scalar, Wrapper> column(size_t const index) const;
      Vector<Scalar, Wrapper> solve(Vector<Scalar, Wrapper> const & rhs);


      DenseMatrix adjoint() const;
      DenseMatrix inverse() const;
      DenseMatrix permute_columns(std::vector<std::size_t> const & permutation) const;
      DenseMatrix operator+(DenseMatrix const & matrix) const;

      typename eigen::DenseBlock<Scalar> block(size_t startRow, size_t startCol, size_t blockRows, size_t blockCols)
      { return wrapped_.block(startRow, startCol, blockRows, blockCols); }

      typename eigen::DenseBlock_const<Scalar> block(size_t startRow, size_t startCol, size_t blockRows, size_t blockCols) const
      { return wrapped_.block(startRow, startCol, blockRows, blockCols); }
      /*typename Wrapper::template DenseBlock<Scalar> block(std::size_t const startRow, 
                                                          std::size_t const startCol, 
                                                          std::size_t const blockRows, 
                                                          std::size_t const blockCols);
      
      typename Wrapper::template DenseBlock_const<Scalar> block(std::size_t const startRow, 
                                                                std::size_t const startCol, 
                                                                std::size_t const blockRows, 
                                                                std::size_t const blockCols) const;
*/

    public:
      typename Wrapper::template DenseMatrix<Scalar> wrapped_;
  };

  //! Null matrix
  template<typename Scalar, typename Wrapper> inline
  DenseMatrix<Scalar, Wrapper> DenseMatrix<Scalar, Wrapper>::Zero(std::size_t const numberOfRows, 
                                                                  std::size_t const numberOfColumns)
  { return DenseMatrix<Scalar>(Wrapper::template Zero<Scalar>(numberOfRows, numberOfColumns)); }

  //! Identity matrix
  template<typename Scalar, typename Wrapper> inline
  DenseMatrix<Scalar, Wrapper> DenseMatrix<Scalar, Wrapper>::Identity(std::size_t const dimension)
  { return DenseMatrix<Scalar, Wrapper>(Wrapper::template Identity<Scalar>(dimension)); }

  // TODO: implement a version which does not depend on eigen
  //! Construction from a sparse matrix
  template<typename Scalar, typename Wrapper> inline
  DenseMatrix<Scalar, Wrapper>::DenseMatrix(SparseMatrix<Scalar, Wrapper> const & matrix)
  : wrapped_(matrix.wrapped_.template triangularView<Eigen::Upper>())
  { wrapped_ += matrix.wrapped_.transpose().template triangularView<Eigen::StrictlyLower>(); }

  //! Construction from a matrix block
  template<typename Scalar, typename Wrapper> inline
  DenseMatrix<Scalar, Wrapper>::DenseMatrix(typename Wrapper::template DenseBlock_const<Scalar> const & block)
  : wrapped_(block)
  {}

  // TODO: move this function to Vector class
  //! Construction from a vector
  template<typename Scalar, typename Wrapper> inline
  DenseMatrix<Scalar, Wrapper>::DenseMatrix(Vector<Scalar, Wrapper> u, std::size_t const numberOfRows, std::size_t const numberOfColumns)
  : DenseMatrix<Scalar, Wrapper>(numberOfRows, numberOfColumns)
  {
    if (u.size() != numberOfColumns * numberOfRows)
      throw std::invalid_argument("Matrix reshape : sizes incoherent !");

    for (std::size_t j = 0; j < numberOfColumns; j++)
      for (std::size_t i = 0; i < numberOfRows; i++)
        wrapped_(i, j) = u(i * u.size() + j);
  }

  //! Construction with prescribed size
  template<typename Scalar, typename Wrapper> inline
  DenseMatrix<Scalar, Wrapper>::DenseMatrix(std::size_t numberOfRows, std::size_t numberOfColumns)
  : wrapped_(numberOfRows, numberOfColumns)
  {}

  //! Construction from wrapped DenseMatrix
  template<typename Scalar, typename Wrapper> inline
  DenseMatrix<Scalar, Wrapper>::DenseMatrix(typename Wrapper::template DenseMatrix<Scalar> const & matrix)
  : wrapped_(matrix)
  {}

  //! Fill with a scalar
  template<typename Scalar, typename Wrapper> inline
  void DenseMatrix<Scalar, Wrapper>::fill(Scalar const scalar)
  { wrapped_.fill(scalar); }
 
  //! Returns the number of rows
  template<typename Scalar, typename Wrapper> inline
  std::size_t DenseMatrix<Scalar, Wrapper>::number_of_rows() const
  { return Wrapper::number_of_rows(wrapped_); }

  //! Returns the number of columns
  template<typename Scalar, typename Wrapper> inline
  std::size_t DenseMatrix<Scalar, Wrapper>::number_of_columns() const
  { return Wrapper::number_of_columns(wrapped_); }
  
  //! Returns the determinant
  template<typename Scalar, typename Wrapper> inline
  Scalar DenseMatrix<Scalar, Wrapper>::determinant() const
  { return Wrapper::determinant(wrapped_); }

  // TODO: unstable version -> should probably return the ratio of absolutes values (but used like this in hartree_fock)
  //! Returns reciprocal condition number
  template<typename Scalar, typename Wrapper> inline
  Scalar DenseMatrix<Scalar, Wrapper>::rcond() const
  {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(wrapped_);

    Vector<Scalar, Wrapper> D = svd.singularValues();
    Scalar lmin = D.min();
    Scalar lmax = D.max();

    return lmin / lmax;
  }

  //! Returns the trace
  template<typename Scalar, typename Wrapper> inline
  Scalar DenseMatrix<Scalar, Wrapper>::trace() const
  { return Wrapper::trace(wrapped_); }

  //! Returns a mutable coefficient
  template<typename Scalar, typename Wrapper> inline
  Scalar & DenseMatrix<Scalar, Wrapper>::operator()(size_t const rowIndex, size_t const columnIndex)
  { return wrapped_(rowIndex, columnIndex); }

  //! Returns an immutable coefficient
  template<typename Scalar, typename Wrapper> inline
  Scalar const & DenseMatrix<Scalar, Wrapper>::operator()(size_t const rowIndex, size_t const columnIndex) const
  { return wrapped_(rowIndex, columnIndex); }

  //! Returns a column
  template<typename Scalar, typename Wrapper> inline
  Vector<Scalar, Wrapper> DenseMatrix<Scalar, Wrapper>::column(size_t const index) const
  { return Vector<Scalar, Wrapper>(wrapped_.col(index)); }

  //! Solve a linear system
  template<typename Scalar, typename Wrapper> inline
  Vector<Scalar, Wrapper> DenseMatrix<Scalar, Wrapper>::solve(Vector<Scalar, Wrapper> const & rhs)
  {
    Vector<Scalar> sol(number_of_rows());
    sol.wrapped_ = wrapped_.partialPivLu().solve(rhs.wrapped_);
    return sol;
  }


  //! Returns the adjoint matrix
  template<typename Scalar, typename Wrapper> inline
  DenseMatrix<Scalar, Wrapper> DenseMatrix<Scalar, Wrapper>::adjoint() const
  { return DenseMatrix<Scalar>(Wrapper::adjoint(wrapped_)); }

  // TODO: write a non-pessimized version
  // with CwiseBinaryOp from Eigen
  //! Sum of matrices
  template<typename Scalar, typename Wrapper> inline
  DenseMatrix<Scalar, Wrapper> DenseMatrix<Scalar, Wrapper>::operator+(DenseMatrix const & other) const
  { return typename Wrapper::template DenseMatrix<Scalar>(this->wrapped_ + other.wrapped_); }

  //! Permute columns
  template<typename Scalar, typename Wrapper> inline
  DenseMatrix<Scalar, Wrapper> DenseMatrix<Scalar, Wrapper>::permute_columns(std::vector<std::size_t> const & permutation) const
  {
    DenseMatrix<Scalar, Wrapper> permuted(number_of_rows(), permutation.size());
    for (size_t i = 0; i < permutation.size(); ++i)
      permuted.wrapped_.col(i) = wrapped_.col(permutation[i]);
    return permuted;
  }

  //! Returns the inverse
  template<typename Scalar, typename Wrapper> inline
  DenseMatrix<Scalar, Wrapper> DenseMatrix<Scalar, Wrapper>::inverse() const
  { return DenseMatrix(Wrapper::inverse(wrapped_)); }

  //! Returns a mutable submatrix
  /*template<typename Scalar, typename Wrapper> inline
  typename Wrapper::template DenseBlock<Scalar> DenseMatrix<Scalar, Wrapper>::block(std::size_t const startRow, 
                                                          std::size_t const startCol, 
                                                          std::size_t const blockRows, 
                                                          std::size_t const blockCols)
  { return Wrapper::block(wrapped_, startRow, startCol, blockRows, blockCols); }

  //! Returns an immutable submatrix
  template<typename Scalar, typename Wrapper> inline
  typename Wrapper::template DenseBlock_const<Scalar> DenseMatrix<Scalar, Wrapper>::block(std::size_t const startRow, 
                                                          std::size_t const startCol, 
                                                          std::size_t const blockRows, 
                                                          std::size_t const blockCols) const
  { return Wrapper::block(wrapped_, startRow, startCol, blockRows, blockCols); }*/

  //! Matrix-vector product
  template<typename Scalar, typename Wrapper> inline
  Vector<Scalar, Wrapper> operator*(DenseMatrix<Scalar, Wrapper> const & matrix, 
                                    Vector<Scalar, Wrapper> const & vector)
  { return Vector<Scalar, Wrapper>(matrix.wrapped_ * vector.wrapped_); }
  
  //! Matrix-matrix product
  template<typename Scalar, typename Wrapper> inline
  DenseMatrix<Scalar, Wrapper> operator*(DenseMatrix<Scalar, Wrapper> const & lhs, 
                                         DenseMatrix<Scalar, Wrapper> const & rhs)
  { return DenseMatrix<Scalar, Wrapper>(lhs.wrapped_ * rhs.wrapped_); }

  //! Piecewise division
  template<typename Scalar, typename Wrapper> inline
  DenseMatrix<Scalar, Wrapper> piecewiseDivision(DenseMatrix<Scalar, Wrapper> const & lhs,
                                                 DenseMatrix<int, Wrapper> const & rhs)
  {
    if (lhs.number_of_rows() != rhs.number_of_rows() || 
        lhs.number_of_columns() != rhs.number_of_columns())
      throw std::invalid_argument("Can only divide matrices of same size !");
    
    DenseMatrix<Scalar, Wrapper> division(lhs.number_of_rows(), rhs.number_of_columns());
    for (std::size_t j = 0; j < lhs.number_of_columns(); j++)
      for (std::size_t i = 0; i < lhs.number_of_rows(); i++)
        division(i, j) = lhs(i, j) / rhs(i, j);
    return division;
  }

}

#endif
