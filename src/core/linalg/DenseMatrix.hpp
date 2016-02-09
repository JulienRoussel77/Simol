#ifndef SIMOL_DENSEMATRIX_HPP
#define SIMOL_DENSEMATRIX_HPP

#include "core/io/MatrixMarketFile.hpp"

#include <Eigen/Dense>

#include "eigen.hpp"

#include "Vector.hpp"

namespace simol
{

  template<class ScalarType, template<class> class WrappedLibrary = eigen>
  class DenseMatrix;

  template<class ScalarType>
  class DenseMatrix<ScalarType,eigen>
  {
      public:

          DenseMatrix(DenseMatrix const & matrix) = default;

          DenseMatrix(std::size_t numberOfRows,
                      std::size_t numberOfColumns);

          DenseMatrix(typename eigen<ScalarType>::DenseMatrixType const & wrappedMatrix);

          Vector<ScalarType>
          column(size_t const index);

          DenseMatrix(MatrixMarketFile const & file);

          ScalarType const &
          operator()(size_t const rowIndex,
                     size_t const columnIndex) const;

          ScalarType &
          operator()(size_t const rowIndex,
                     size_t const columnIndex);

          size_t
          numberOfRows() const;

          size_t
          numberOfColumns() const;

          double
          rcond() const;

      public:
          typename eigen<ScalarType>::DenseMatrixType wrapped_;
  };

  template<class ScalarType>
  DenseMatrix<ScalarType, eigen>::DenseMatrix(std::size_t numberOfRows,
              std::size_t numberOfColumns)
  : wrapped_(numberOfRows, numberOfColumns)
  {}

  template<class ScalarType>
  DenseMatrix<ScalarType, eigen>::DenseMatrix(typename eigen<ScalarType>::DenseMatrixType const & wrappedMatrix)
  : wrapped_(wrappedMatrix)
  {}

  template<class ScalarType> Vector<ScalarType>
  DenseMatrix<ScalarType, eigen>::column(size_t const index)
  {
      Vector<ScalarType> columnVector(wrapped_.rows());
      columnVector.wrapped_.col(index);
      return columnVector;
  }

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

  template<class ScalarType> ScalarType const &
  DenseMatrix<ScalarType, eigen>::operator()(size_t const rowIndex,
                                             size_t const columnIndex) const
  { return wrapped_(rowIndex, columnIndex); }

  template<class ScalarType> ScalarType &
  DenseMatrix<ScalarType, eigen>::operator()(size_t const rowIndex,
                                             size_t const columnIndex)
  { return wrapped_(rowIndex, columnIndex); }

  template<class ScalarType> size_t
  DenseMatrix<ScalarType, eigen>::numberOfRows() const
  { return wrapped_.rows(); }

  template<class ScalarType> size_t
  DenseMatrix<ScalarType, eigen>::numberOfColumns() const
  { return wrapped_.cols(); }

  template<class ScalarType> double
  DenseMatrix<ScalarType, eigen>::rcond() const
  {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(wrapped_);

    Vector<double> D = svd.singularValues();
    double lmin = D.min();
    double lmax = D.max();

    return lmin/lmax;
  }

}

#endif
