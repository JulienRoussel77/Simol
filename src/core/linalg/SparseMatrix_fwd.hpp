#ifndef SIMOL_SPARSEMATRIX_FWD_HPP
#define SIMOL_SPARSEMATRIX_FWD_HPP

#include "eigen.hpp"

namespace simol
{
  template<class ScalarType, template<class> class Library = eigen>
  class SparseMatrix;
}

#endif
