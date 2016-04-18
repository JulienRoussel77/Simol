#ifndef EIGEN_HPP
#define EIGEN_HPP

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wsign-compare"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCore>
#include <Eigen/SVD>
#pragma GCC diagnostic pop


namespace simol
{

  template<class ScalarType>
  struct eigen
  {
    typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> VectorType;
    typedef Eigen::SparseMatrix<ScalarType> SparseMatrixType;
    typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> DenseMatrixType;
    typedef Eigen::Map<VectorType> VectorMap;
    typedef Eigen::JacobiSVD<DenseMatrixType> SVDType;
    typedef typename DenseMatrixType::AdjointReturnType AdjointReturnType;
    typedef typename DenseMatrixType::ColXpr DenseMatrixColumn;
    typedef Eigen::Block<DenseMatrixType> DenseBlock;
    typedef Eigen::Block<DenseMatrixType const> DenseBlock_const;

    static std::size_t length(VectorType const & vector)
    {
      return 0;
    }
  };

}

#endif
