#ifndef SIMOL_DENSEMATRIX_HPP
#define SIMOL_DENSEMATRIX_HPP

#include "simol/core/io/MatrixMarketFile.hpp"
#include "simol/core/linalg/Vector.hpp"

#include <ostream>

namespace simol
{

  template<class ScalarType, template<class> class WrappedLibrary = eigen>
  class DenseMatrix;

}

#include "simol/core/linalg/DenseMatrix_eigen.hpp"

#endif
