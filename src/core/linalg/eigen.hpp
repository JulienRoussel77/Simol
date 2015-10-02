#ifndef EIGEN_HPP
#define EIGEN_HPP

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wsign-compare"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#pragma GCC diagnostic pop

namespace simol
{

  template<class ScalarType>
  struct eigen
  {
    typedef Eigen::Matrix<ScalarType,Eigen::Dynamic,1> VectorType;
    typedef Eigen::SparseMatrix<ScalarType> SparseMatrixType;
    typedef Eigen::Map<VectorType> VectorMap;
  };

}

#endif
