#include "Vector.hpp"

namespace simol {
    
  Vector<double,eigen> operator*(double const& lambda, Vector<double,eigen> const& v)
    {
      return v*lambda;
    }
    
    double dot(Vector<double,eigen> const& u, Vector<double,eigen> const& v)
    {
      return u.wrapped_.dot(v.wrapped_);
    }
}


  