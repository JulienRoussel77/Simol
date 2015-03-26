#ifndef EIGEN_HPP
#define EIGEN_HPP

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Core>

namespace simol
{

  template<class ScalarType>
  struct eigen
  {
    typedef Eigen::Matrix<ScalarType,Eigen::Dynamic,1> VectorType;
    typedef Eigen::SparseMatrix<ScalarType> SparseMatrixType;
  };

}

#endif