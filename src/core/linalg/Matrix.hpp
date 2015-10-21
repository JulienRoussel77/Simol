#ifndef SIMOL_MATRIX_HPP
#define SIMOL_MATRIX_HPP

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
          
          DenseMatrix(typename eigen<ScalarType>::DenseMatrixType const & wrappedMatrix)
          : wrapped_(wrappedMatrix)
          {}
          
          //Vector<ScalarType> column(size_t const index)
          //{ return wrapped_.col(index); }
          
          /*typename eigen<ScalarType>::AdjointReturnType 
          adjoint() const
          { return wrapped_.adjoint(); }*/
          
          ScalarType const & operator()(size_t const rowIndex,
                                        size_t const columnIndex) const
          { return wrapped_(rowIndex, columnIndex); }
          
          ScalarType & operator()(size_t const rowIndex,
                                  size_t const columnIndex)
          { return wrapped_(rowIndex, columnIndex); }
          
          size_t number_of_rows() const
          { return wrapped_.rows(); }
          
          size_t number_of_columns() const
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
