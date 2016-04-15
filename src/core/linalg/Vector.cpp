#include "simol/core/linalg/Vector.hpp"

namespace simol {

  Vector<double,eigen> operator*(double const& lambda, Vector<double,eigen> const& v)
    {
      return v*lambda;
    }

    double dot(Vector<double,eigen> const& u, Vector<double,eigen> const& v)
    {
      return u.wrapped_.dot(v.wrapped_);
    }

    /*Vector<double,eigen> piecewiseDivision(Vector<double,eigen> const& u, Vector<double,eigen> const& v)
    {
      if (u.size() != v.size())
        throw std::invalid_argument("Can only divide vectors of same size !");
      Vector<double,eigen> w(u.size());
      for (int i=0; i<w.size(); i++)
        w(i) = u(i) / v(i);
      return w;
    }*/

  Vector<double,eigen> piecewiseDivision(Vector<double,eigen> const& u, Vector<double,eigen> const& v)
  {
    if (u.size() != v.size())
      throw std::invalid_argument("Can only divide vectors of same size !");
    Vector<double,eigen> w(u.size());
    for (size_t i=0; i<u.size(); i++)
        w(i) = u(i) / v(i);
    return w;
  }

  Vector<double,eigen> piecewiseDivision(Vector<double,eigen> const& u, Vector<size_t,eigen> const& v)
  {
    if (u.size() != v.size())
      throw std::invalid_argument("Can only divide vectors of same size !");
    Vector<double,eigen> w(u.size());
    for (size_t i=0; i<u.size(); i++)
        w(i) = u(i) / v(i);
    return w;
  }
}



