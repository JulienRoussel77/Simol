#include "Vector.hpp"

namespace simol {
    
  Vector<double,eigen> operator*(double const& lambda, Vector<double,eigen> const& v)
    {
      return v*lambda;
    }
}


  