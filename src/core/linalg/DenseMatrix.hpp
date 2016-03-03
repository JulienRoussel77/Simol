#ifndef SIMOL_DENSEMATRIX_HPP
#define SIMOL_DENSEMATRIX_HPP

#include "core/io/MatrixMarketFile.hpp"


#include "Vector.hpp"

#include <ostream>

namespace simol
{

  template<class ScalarType, template<class> class WrappedLibrary = eigen>
  class DenseMatrix;

}

#include "DenseMatrix_eigen.hpp"

#endif
