#ifndef SIMOL_VECTOR_HPP
#define SIMOL_VECTOR_HPP

#include "VectorInterface.hpp"

#include "eigen.hpp"

namespace simol
{
  template<class ScalarType = double, template<class> class WrappedLibrary = eigen>
  class Vector;
}

#include "Vector_eigen.hpp"

#endif
