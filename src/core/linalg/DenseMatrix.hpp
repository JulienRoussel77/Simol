#ifndef SIMOL_DENSEMATRIX_HPP
#define SIMOL_DENSEMATRIX_HPP

#include "MatrixMarketFile.hpp"

#include <Eigen/Dense>

#include "eigen.hpp"

#include "Vector.hpp"

namespace simol
{
  //template<class ScalarType>
  //using DenseMatrix = Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>;
  
  template<class ScalarType, template<class> class WrappedLibrary = eigen>
  class DenseMatrix;
  
  template<class ScalarType>
  class DenseMatrix<ScalarType,eigen>
  {
      public:
          DenseMatrix(std::size_t numberOfRows,
                      std::size_t numberOfColumns)
          : wrapped_(numberOfRows, numberOfColumns)
          {}
          
          DenseMatrix(DenseMatrix const & matrix) = default;
          
          DenseMatrix(typename eigen<ScalarType>::DenseMatrixType const & wrappedMatrix)
          : wrapped_(wrappedMatrix)
          {}
          
          Vector<ScalarType> column(size_t const index)
          { 
              Vector<ScalarType> columnVector(wrapped_.rows());
              columnVector.wrapped_.col(index); 
              return columnVector;
          }
          
          
          DenseMatrix(MatrixMarketFile const & file)
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
          
          /*typename eigen<ScalarType>::AdjointReturnType 
          adjoint() const
          { return wrapped_.adjoint(); }*/
          
          ScalarType const & operator()(size_t const rowIndex,
                                        size_t const columnIndex) const
          { return wrapped_(rowIndex, columnIndex); }
          
          ScalarType & operator()(size_t const rowIndex,
                                  size_t const columnIndex)
          { return wrapped_(rowIndex, columnIndex); }
          
          size_t numberOfRows() const
          { return wrapped_.rows(); }
          
          size_t numberOfColumns() const
          { return wrapped_.cols(); }
          
          double rcond() const
          {
            Eigen::JacobiSVD<Eigen::MatrixXd> svd(wrapped_);

            Vector<double> D = svd.singularValues();
            double lmin = D.min();
            double lmax = D.max();

            return lmin/lmax;
          }
          
      public:
          typename eigen<ScalarType>::DenseMatrixType wrapped_;
  };
  

}

#endif
