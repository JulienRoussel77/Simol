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

      void fill(Scalar scalar);

      double rcond() const;


      Vector<Scalar> column(size_t const index) const;






      // TODO: write a non-pessimized version (return type may different in eigen)
      DenseMatrix adjoint() const
      { return DenseMatrix(wrapped_.adjoint()); }

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

  template<class Scalar> inline
  DenseMatrix<Scalar, eigen> DenseMatrix<Scalar, eigen>::operator*(Scalar const& scalar) const
  {
    DenseMatrix<Scalar> u(*this);
    u.wrapped_ *= scalar;
    return u;
  }

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


  template<class Scalar>
  DenseMatrix<Scalar, eigen> operator*(Scalar const& scalar, DenseMatrix<Scalar, eigen> const& A)
  {return A * scalar;}


}

#endif

