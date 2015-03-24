#ifndef VECTORWRAPPER_HPP
#define VECTORWRAPPER_HPP

#include "stl.hpp"
#include "eigen.hpp"

#include "VectorWrapper_fwd.hpp"

namespace simol
{

  template<class ScalarType>
  class VectorWrapper<ScalarType,stl> : public stl<ScalarType>
  {};

  template<class ScalarType>
  class VectorWrapper<ScalarType,eigen> : public eigen<ScalarType>
  {};

}

#endif
