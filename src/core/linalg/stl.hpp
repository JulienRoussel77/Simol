#ifndef SIMOL_STL_HPP
#define SIMOL_STL_HPP

#include <vector>

namespace simol
{

  template<class ScalarType>
  struct stl
  {
    typedef std::vector<ScalarType> VectorType;
  };

}

#endif
