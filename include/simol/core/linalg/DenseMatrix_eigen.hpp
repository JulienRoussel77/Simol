#ifndef SIMOL_DENSEMATRIX_EIGEN_HPP
#define SIMOL_DENSEMATRIX_EIGEN_HPP

#include <Eigen/Core>
#include <Eigen/Dense>

#include "eigen.hpp"
#include "SparseMatrix.hpp"

namespace simol
{
  template<class Scalar>
  class DenseMatrix<Scalar, eigen>
  {

    public:

      static DenseMatrix Zero(std::size_t const numberOfRows, std::size_t numberOfColumns);
      static DenseMatrix Identity(std::size_t const dimension);

      DenseMatrix(DenseMatrix const & matrix) = default;
      DenseMatrix(MatrixMarketFile const & file);
      DenseMatrix(SparseMatrix<Scalar, eigen> const & matrix);
      DenseMatrix(std::size_t numberOfRows, std::size_t numberOfColumns);
      DenseMatrix(typename eigen::DenseMatrix<Scalar> const & wrappedMatrix);
      DenseMatrix(typename eigen::DenseBlock_const<Scalar> const & block);
      DenseMatrix(Vector<double, eigen> u, std::size_t numberOfRows, std::size_t numberOfColumns);

      void fill(Scalar scalar);

      std::size_t number_of_rows() const;
      std::size_t number_of_columns() const;

      double rcond() const;

      Scalar trace() const;
      Scalar determinant() const;

      DenseMatrix inverse() const;

      Vector<Scalar> column(size_t const index) const;
      Vector<Scalar> solve(Vector<Scalar> const & rhs)
      {
        Vector<Scalar> sol(number_of_rows());
        sol.wrapped_ = wrapped_.partialPivLu().solve(rhs.wrapped_);
        return sol;
      }


      typename eigen::DenseBlock<Scalar> block(size_t startRow, size_t startCol, size_t blockRows, size_t blockCols)
      { return wrapped_.block(startRow, startCol, blockRows, blockCols); }

      typename eigen::DenseBlock_const<Scalar> block(size_t startRow, size_t startCol, size_t blockRows, size_t blockCols) const
      { return wrapped_.block(startRow, startCol, blockRows, blockCols); }

      DenseMatrix permute_columns(std::vector<std::size_t> const & permutation) const
      {
        DenseMatrix permuted(number_of_rows(), permutation.size());
        for (size_t i = 0; i < permutation.size(); ++i)
          permuted.wrapped_.col(i) = wrapped_.col(permutation[i]);
        return permuted;
      }


      Scalar const & operator()(size_t const rowIndex, size_t const columnIndex) const;
      Scalar & operator()(size_t const rowIndex, size_t const columnIndex);




      // TODO: write a non-pessimized version
      // with CwiseBinaryOp from Eigen
      DenseMatrix operator+(DenseMatrix const & other) const
      { return WrappedType(this->wrapped_ + other.wrapped_); }

      // TODO: write a non-pessimized version (return type may different in eigen)
      DenseMatrix adjoint() const
      { return DenseMatrix(wrapped_.adjoint()); }

      DenseMatrix operator*(DenseMatrix const & matrix) const
      { return DenseMatrix( WrappedType(wrapped_ * matrix.wrapped_) ); }

      DenseMatrix operator*=(Scalar const& scalar)
      {
        wrapped_ *= scalar;
        return *this;
      }

      DenseMatrix<Scalar, eigen> operator*(Scalar const& scalar) const;


      DenseMatrix & operator*=(Scalar const scalar)
      {
        wrapped_ *= scalar;
        return *this;
      }

    public:
      typedef typename eigen::DenseMatrix<Scalar> WrappedType;
      typename eigen::DenseMatrix<Scalar> wrapped_;
  };

  template<typename Scalar> inline
  DenseMatrix<Scalar, eigen> DenseMatrix<Scalar, eigen>::Zero(std::size_t numberOfRows, size_t numberOfColumns)
  { return DenseMatrix(eigen::Zero<Scalar>(numberOfRows, numberOfColumns)); }

  template<typename Scalar> inline
  DenseMatrix<Scalar, eigen> DenseMatrix<Scalar, eigen>::Identity(std::size_t const dimension)
  { return DenseMatrix(eigen::Identity<Scalar>(dimension)); }

  //! Construction from a sparse matrix
  template<typename Scalar> inline
  DenseMatrix<Scalar, eigen>::DenseMatrix(SparseMatrix<Scalar, eigen> const & matrix)
    : wrapped_(matrix.wrapped_.template triangularView<Eigen::Upper>())
  { wrapped_ += matrix.wrapped_.transpose().template triangularView<Eigen::StrictlyLower>(); }

  //! Construction from a matrix block
  template<typename Scalar> inline
  DenseMatrix<Scalar, eigen>::DenseMatrix(typename eigen::DenseBlock_const<Scalar> const & block)
    : wrapped_(block)
  {}

  //! Construction with prescribed size
  template<class Scalar>
  DenseMatrix<Scalar, eigen>::DenseMatrix(std::size_t numberOfRows, std::size_t numberOfColumns)
    : wrapped_(numberOfRows, numberOfColumns)
  {}

  //! Construction from eigen's DenseMatrix
  template<class Scalar>
  DenseMatrix<Scalar, eigen>::DenseMatrix(typename eigen::DenseMatrix<Scalar> const & wrappedMatrix)
    : wrapped_(wrappedMatrix)
  {}

  //! Returns the trace
  template<typename Scalar> inline
  Scalar DenseMatrix<Scalar, eigen>::trace() const
  { return eigen::trace(wrapped_); }

  //! Returns the determinant
  template<typename Scalar> inline
  Scalar DenseMatrix<Scalar, eigen>::determinant() const
  { return eigen::determinant(wrapped_); }

  //! Returns the inverse
  template<typename Scalar> inline
  DenseMatrix<Scalar, eigen> DenseMatrix<Scalar, eigen>::inverse() const
  { return DenseMatrix(eigen::inverse(wrapped_)); }

  //! Returns a column
  template<class Scalar> Vector<Scalar>
  DenseMatrix<Scalar, eigen>::column(size_t const index) const
  { return Vector<Scalar>(wrapped_.col(index)); }


  //! Construction from a Matrix Market file
  template<class Scalar>
  DenseMatrix<Scalar, eigen>::DenseMatrix(MatrixMarketFile const & file)
    : wrapped_(eigen::Zero<Scalar>(file.numberOfRows(), file.numberOfRows()))
  {
    std::vector< Eigen::Triplet<Scalar, std::size_t> > nonzeros(file.numberOfNonzeros());
    for(size_t nonzeroIndex = 0; nonzeroIndex < file.numberOfNonzeros(); ++nonzeroIndex)
    {
      int rowIndex;
      int columnIndex;
      Scalar nonzero;
      fscanf(file.content(), "%d %d %lg\n", &rowIndex, &columnIndex, &nonzero);
      wrapped_(rowIndex, columnIndex) = nonzero;
    }
  }

  // TODO: move this function to Vector class
  //! Construction from a vector
  template<class Scalar>
  DenseMatrix<Scalar, eigen>::DenseMatrix(Vector<double, eigen> u, std::size_t numberOfRows, std::size_t numberOfColumns):
    DenseMatrix(numberOfRows, numberOfColumns)
  {
    if (u.size() != numberOfColumns * numberOfRows)
      throw std::invalid_argument("Matrix reshape : sizes incoherent !");

    for (std::size_t j = 0; j < numberOfColumns; j++)
      for (std::size_t i = 0; i < numberOfRows; i++)
        wrapped_(i, j) = u(i * u.size() + j);
  }

  //! Returns a coefficient
  template<class Scalar> inline
  Scalar const & DenseMatrix<Scalar, eigen>::operator()(size_t const rowIndex, size_t const columnIndex) const
  { return wrapped_(rowIndex, columnIndex); }

  //! Modify a coefficient
  template<class Scalar> inline
  Scalar & DenseMatrix<Scalar, eigen>::operator()(size_t const rowIndex, size_t const columnIndex)
  { return wrapped_(rowIndex, columnIndex); }

  //! Returns the number of rows
  template<typename Scalar> inline
  std::size_t DenseMatrix<Scalar, eigen>::number_of_rows() const
  { return eigen::number_of_rows(wrapped_); }

  //! Returns the number of columns
  template<typename Scalar> inline
  std::size_t DenseMatrix<Scalar, eigen>::number_of_columns() const
  { return eigen::number_of_rows(wrapped_); }

  template<class Scalar> inline
  DenseMatrix<Scalar, eigen> DenseMatrix<Scalar, eigen>::operator*(Scalar const& scalar) const
  {
    DenseMatrix<Scalar> u(*this);
    u.wrapped_ *= scalar;
    return u;
  }

  //! Returns reciprocal condition number
  template<class Scalar> inline
  double DenseMatrix<Scalar, eigen>::rcond() const
  {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(wrapped_);

    Vector<double> D = svd.singularValues();
    double lmin = D.min();
    double lmax = D.max();

    return lmin / lmax;
  }

  template<class Scalar> inline
  void DenseMatrix<Scalar, eigen>::fill(Scalar scalar)
  { wrapped_.fill(scalar); }

  //! Print a matrix into a file
  template<class Scalar> inline
  std::ostream & operator<<(std::ostream & output, DenseMatrix<Scalar, eigen> const & matrixToPrint)
  { return output << matrixToPrint.wrapped_; }

  //! Returns matrix-vector product
  template<typename Scalar> inline
  Vector<Scalar> operator*(DenseMatrix<Scalar, eigen> const & matrix, Vector<Scalar> const & vector)
  { return Vector<Scalar>(matrix.wrapped_ * vector.wrapped_); }

  // TODO: rename it Identity and call the right function from eigen
  template<class Scalar>
  DenseMatrix<Scalar, eigen> eye(size_t const nbOfRows, size_t const nbOfColumns)
  {
    DenseMatrix<Scalar, eigen> A(nbOfRows, nbOfColumns);
    for (std::size_t i = 0; i < std::min(nbOfRows, nbOfColumns); i++)
      A(i, i) = 1;
    //A.wrapped_.setIdentity();
    return A;
  }

  template<class ScalarA, class ScalarB>
  DenseMatrix<double, eigen> piecewiseDivision(DenseMatrix<ScalarA, eigen> const& A, DenseMatrix<ScalarB, eigen> const& B)
  {
    if (A.number_of_rows() != B.number_of_rows() || A.number_of_columns() != B.number_of_columns())
      throw std::invalid_argument("Can only divide matrices of same size !");
    DenseMatrix<double, eigen> C(A.number_of_rows(), B.number_of_columns());
    for (std::size_t j = 0; j < A.number_of_columns(); j++)
      for (std::size_t i = 0; i < A.number_of_rows(); i++)
        C(i, j) = A(i, j) / B(i, j);
    return C;
  }


  template<class Scalar>
  DenseMatrix<Scalar, eigen> operator*(Scalar const& scalar, DenseMatrix<Scalar, eigen> const& A)
  {return A * scalar;}


}

#endif

