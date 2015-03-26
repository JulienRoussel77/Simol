#ifndef FROMEIGEN_HPP
#define FROMEIGEN_HPP

#include<Eigen/Sparse>
#include<Eigen/SparseCore>

namespace simol
{

  template<class ScalarType>
  struct fromEigen
  {
    typedef Eigen::SparseMatrix<ScalarType> WrappedType;
  };

}

#endif
