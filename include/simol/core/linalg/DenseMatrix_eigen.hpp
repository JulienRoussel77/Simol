#ifndef SIMOL_DENSEMATRIX_EIGEN_HPP
#define SIMOL_DENSEMATRIX_EIGEN_HPP

#include <Eigen/Core>
#include <Eigen/Dense>

#include "eigen.hpp"
#include "SparseMatrix.hpp"

namespace simol
{
  template<class ScalarType>
  class DenseMatrix<ScalarType, eigen>
  {

    public:

      DenseMatrix(DenseMatrix const & matrix) = default;
      DenseMatrix(SparseMatrix<ScalarType, eigen> const & matrix);
      DenseMatrix(std::size_t numberOfRows, std::size_t numberOfColumns);
      DenseMatrix(typename eigen<ScalarType>::DenseMatrixType const & wrappedMatrix);
      DenseMatrix(typename eigen<ScalarType>::DenseBlock_const const & block);
      DenseMatrix(Vector<double, eigen> u, std::size_t numberOfRows, std::size_t numberOfColumns);

      Vector<ScalarType> solve(Vector<ScalarType> const & rhs)
      {
        Vector<ScalarType> sol(numberOfRows());
        sol.wrapped_ = wrapped_.partialPivLu().solve(rhs.wrapped_);
        return sol;
      }

      Vector<ScalarType>
      column(size_t const index) const;

      typename eigen<ScalarType>::DenseBlock block(size_t startRow, size_t startCol, size_t blockRows, size_t blockCols)
      { return wrapped_.block(startRow, startCol, blockRows, blockCols); }

      typename eigen<ScalarType>::DenseBlock_const block(size_t startRow, size_t startCol, size_t blockRows, size_t blockCols) const
      { return wrapped_.block(startRow, startCol, blockRows, blockCols); }

      DenseMatrix permute_columns(std::vector<std::size_t> const & permutation) const
      {
        DenseMatrix permuted(numberOfRows(), permutation.size());
        for (size_t i = 0; i < permutation.size(); ++i)
          permuted.wrapped_.col(i) = wrapped_.col(permutation[i]);
        return permuted;
      }

      DenseMatrix(MatrixMarketFile const & file);

      ScalarType const & operator()(size_t const rowIndex, size_t const columnIndex) const;
      ScalarType & operator()(size_t const rowIndex, size_t const columnIndex);

      size_t numberOfRows() const;
      size_t numberOfColumns() const;

      double rcond() const;

      void fill(ScalarType scalar);

      static DenseMatrix Zero(std::size_t numberOfRows, size_t numberOfColumns)
      { return DenseMatrix(WrappedType::Zero(numberOfRows, numberOfColumns)); }

      static DenseMatrix Identity(std::size_t n)
      { return DenseMatrix(WrappedType::Identity(n)); }

      // TODO: write a non-pessimized version
      // with CwiseBinaryOp from Eigen
      DenseMatrix operator+(DenseMatrix const & other) const
      { return WrappedType(this->wrapped_ + other.wrapped_); }

      // TODO: write a non-pessimized version (return type may different in eigen)
      DenseMatrix adjoint() const
      { return DenseMatrix(wrapped_.adjoint()); }

      DenseMatrix operator*(DenseMatrix const & matrix) const
      { return DenseMatrix( WrappedType(wrapped_ * matrix.wrapped_) ); }

      DenseMatrix operator*=(ScalarType const& scalar)
      {
        wrapped_ *= scalar;
        return *this;
      }

      DenseMatrix<ScalarType, eigen> operator*(ScalarType const& scalar) const;

      ScalarType trace() const
      { return wrapped_.trace(); }

      ScalarType determinant() const
      { return wrapped_.determinant(); }

      DenseMatrix inverse() const
      { return DenseMatrix(WrappedType(wrapped_.inverse())); }


      DenseMatrix & operator*=(ScalarType const scalar)
      {
        wrapped_ *= scalar;
        return *this;
      }

    public:
      typedef typename eigen<ScalarType>::DenseMatrixType WrappedType;
      typename eigen<ScalarType>::DenseMatrixType wrapped_;
  };

  /*//! Construction from a sparse matrix
  template<typename ScalarType> inline
  DenseMatrix<ScalarType, eigen>::DenseMatrix(SparseMatrix<ScalarType, eigen> const & matrix)
    : wrapped_(matrix.wrapped_.template triangularView<Eigen::Upper>())
  { wrapped_ += matrix.wrapped_.transpose().template triangularView<Eigen::StrictlyLower>(); }*/
  
  //! Construction from a sparse matrix
  template<typename ScalarType> inline
  DenseMatrix<ScalarType, eigen>::DenseMatrix(SparseMatrix<ScalarType, eigen> const & matrix)
    : wrapped_(matrix.wrapped_)
  {}

  //! Construction from a matrix block
  template<typename ScalarType> inline
  DenseMatrix<ScalarType, eigen>::DenseMatrix(typename eigen<ScalarType>::DenseBlock_const const & block)
    : wrapped_(block)
  {}

  //! Construction with prescribed size
  template<class ScalarType>
  DenseMatrix<ScalarType, eigen>::DenseMatrix(std::size_t numberOfRows, std::size_t numberOfColumns)
    : wrapped_(numberOfRows, numberOfColumns)
  {}

  //! Construction from eigen's DenseMatrix
  template<class ScalarType>
  DenseMatrix<ScalarType, eigen>::DenseMatrix(typename eigen<ScalarType>::DenseMatrixType const & wrappedMatrix)
    : wrapped_(wrappedMatrix)
  {}

  //! Returns a column
  template<class ScalarType> Vector<ScalarType>
  DenseMatrix<ScalarType, eigen>::column(size_t const index) const
  { return Vector<ScalarType>(wrapped_.col(index)); }


  //! Construction from a Matrix Market file
  template<class ScalarType>
  DenseMatrix<ScalarType, eigen>::DenseMatrix(MatrixMarketFile const & file)
    : wrapped_(eigen<double>::DenseMatrixType::Zero(file.numberOfRows(), file.numberOfColumns()))
  {
    std::vector< Eigen::Triplet<ScalarType, std::size_t> > nonzeros(file.numberOfNonzeros());
    for(size_t nonzeroIndex = 0; nonzeroIndex < file.numberOfNonzeros(); ++nonzeroIndex)
    {
      int rowIndex;
      int columnIndex;
      ScalarType nonzero;
      fscanf(file.content(), "%d %d %lg\n", &rowIndex, &columnIndex, &nonzero);
      wrapped_(rowIndex, columnIndex) = nonzero;
    }
  }

  // TODO: move this function to Vector class
  //! Construction from a vector
  template<class ScalarType>
  DenseMatrix<ScalarType, eigen>::DenseMatrix(Vector<double, eigen> u, std::size_t numberOfRows, std::size_t numberOfColumns):
    DenseMatrix(numberOfRows, numberOfColumns)
  {
    if (u.size() != numberOfColumns * numberOfRows)
      throw std::invalid_argument("Matrix reshape : sizes incoherent !");

    for (std::size_t j = 0; j < numberOfColumns; j++)
      for (std::size_t i = 0; i < numberOfRows; i++)
        wrapped_(i, j) = u(j * numberOfRows+ i);
        //wrapped_(i, j) = u(i * numberOfColumns + j);
  }

  //! Returns a coefficient
  template<class ScalarType> inline
  ScalarType const & DenseMatrix<ScalarType, eigen>::operator()(size_t const rowIndex, size_t const columnIndex) const
  { return wrapped_(rowIndex, columnIndex); }

  //! Modify a coefficient
  template<class ScalarType> inline
  ScalarType & DenseMatrix<ScalarType, eigen>::operator()(size_t const rowIndex, size_t const columnIndex)
  { return wrapped_(rowIndex, columnIndex); }

  //! Returns the number of rows
  template<class ScalarType> inline
  std::size_t DenseMatrix<ScalarType, eigen>::numberOfRows() const
  { return wrapped_.rows(); }

  //! Returns the number of columns
  template<class ScalarType> inline
  std::size_t DenseMatrix<ScalarType, eigen>::numberOfColumns() const
  { return wrapped_.cols(); }

  template<class ScalarType> inline
  DenseMatrix<ScalarType, eigen> DenseMatrix<ScalarType, eigen>::operator*(ScalarType const& scalar) const
  {
    DenseMatrix<ScalarType> u(*this);
    u.wrapped_ *= scalar;
    return u;
  }

  //! Returns reciprocal condition number
  template<class ScalarType> inline
  double DenseMatrix<ScalarType, eigen>::rcond() const
  {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(wrapped_);

    Vector<double> D = svd.singularValues();
    double lmin = D.min();
    double lmax = D.max();

    return lmin / lmax;
  }

  template<class ScalarType> inline
  void DenseMatrix<ScalarType, eigen>::fill(ScalarType scalar)
  { wrapped_.fill(scalar); }

  //! Print a matrix into a file
  template<class ScalarType> inline
  std::ostream & operator<<(std::ostream & output, DenseMatrix<ScalarType, eigen> const & matrixToPrint)
  { return output << matrixToPrint.wrapped_; }

  //! Returns matrix-vector product
  template<typename ScalarType> inline
  Vector<ScalarType> operator*(DenseMatrix<ScalarType, eigen> const & matrix, Vector<ScalarType> const & vector)
  { return Vector<ScalarType>(matrix.wrapped_ * vector.wrapped_); }

  // TODO: rename it Identity and call the right function from eigen
  template<class ScalarType>
  DenseMatrix<ScalarType, eigen> eye(size_t const nbOfRows, size_t const nbOfColumns)
  {
    DenseMatrix<ScalarType, eigen> A(nbOfRows, nbOfColumns);
    for (std::size_t i = 0; i < std::min(nbOfRows, nbOfColumns); i++)
      A(i, i) = 1;
    //A.wrapped_.setIdentity();
    return A;
  }

  // TODO: remove it and call Zero instead
  template<class ScalarType> inline
  DenseMatrix<ScalarType, eigen> zero(size_t const nbOfRows, size_t const nbOfColumns)
  {
    DenseMatrix<ScalarType, eigen> A = DenseMatrix<ScalarType, eigen>::Zero(nbOfRows, nbOfColumns);
    return A;
  }

  template<class ScalarTypeA, class ScalarTypeB>
  DenseMatrix<double, eigen> piecewiseDivision(DenseMatrix<ScalarTypeA, eigen> const& A, DenseMatrix<ScalarTypeB, eigen> const& B)
  {
    if (A.numberOfRows() != B.numberOfRows() || A.numberOfColumns() != B.numberOfColumns())
      throw std::invalid_argument("Can only divide matrices of same size !");
    DenseMatrix<double, eigen> C(A.numberOfRows(), B.numberOfColumns());
    for (std::size_t j = 0; j < A.numberOfColumns(); j++)
      for (std::size_t i = 0; i < A.numberOfRows(); i++)
        C(i, j) = A(i, j) / B(i, j);
    return C;
  }


  template<class ScalarType>
  DenseMatrix<ScalarType, eigen> operator*(ScalarType const& scalar, DenseMatrix<ScalarType, eigen> const& A)
  {return A * scalar;}


}

#endif

