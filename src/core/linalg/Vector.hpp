#ifndef VECTORWRAPPER_HPP
#define VECTORWRAPPER_HPP

#include "eigen.hpp"

namespace simol
{
  template<class ScalarType=double, template<class> class WrappedLibrary=eigen>
  class Vector;

  typedef Vector<double> dvec;

}

#include "Vector_eigen.hpp"

#endif
