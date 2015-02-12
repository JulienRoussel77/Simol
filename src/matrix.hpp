#ifndef SIMOL_MATRIX_HPP
#define SIMOL_MATRIX_HPP

#include <Eigen/Dense>

namespace simol
{

  template<class ScalarType>
  using DenseMatrix = Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>;

  /*template<class ScalarType>
  class DenseMatrix
  {
    public:

    private:
      Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> self_;

  };*/

}

#endif
