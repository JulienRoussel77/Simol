#ifndef VECTORWRAPPER_HPP
#define VECTORWRAPPER_HPP

#include "stl.hpp"
#include "eigen.hpp"
#include <iostream>

namespace simol
{
  template<class ScalarType=double, template<class> class WrappedLibrary=eigen>
  class Vector;
  
  typedef Vector<double> dvec;
  
  //double dot(Vector<double,eigen> const& u, Vector<double,eigen> const& v);
}

namespace simol
{

  template<class ScalarType>
  class Vector<ScalarType,stl> : public stl<ScalarType>
  {};
  

  template<class ScalarType>
  class Vector<ScalarType,eigen>
  {
    //friend double dot(Vector<double,eigen> const& u, Vector<double,eigen> const& v);
    public:
      Vector(size_t const size=0);
      Vector(size_t const size, ScalarType const& lambda);
      Vector(Vector<ScalarType,eigen> const& u);
      size_t size() const;
      ScalarType & operator()(size_t const index);
      ScalarType const & operator()(size_t const index) const;
      ScalarType norm() const;
      Vector<ScalarType,eigen>& fill(ScalarType const& lambda);
      ScalarType dot(Vector<ScalarType,eigen> const& v) const;

      Vector<ScalarType,eigen>& operator+=(Vector<ScalarType,eigen> const& v);
      Vector<ScalarType,eigen>& operator-=(Vector<ScalarType,eigen> const& v);
      Vector<ScalarType,eigen>& operator*=(ScalarType const& lambda);
      Vector<ScalarType,eigen>& operator/=(ScalarType const& lambda);     
      Vector<ScalarType,eigen> operator*(ScalarType const& lambda) const;
      Vector<ScalarType,eigen> operator/(ScalarType const& lambda) const;
      Vector<ScalarType,eigen> operator-() const;
      Vector<ScalarType,eigen> operator+(Vector<ScalarType,eigen> const& v) const;
      Vector<ScalarType,eigen> operator-(Vector<ScalarType,eigen> const& v) const;


    private:
      typename eigen<ScalarType>::VectorType wrapped_;
      
  };
  


  

  Vector<double,eigen> operator*(double const& lambda, Vector<double,eigen> const& v);
  

}

template<class ScalarType, template<class> class WrappedLibrary>
std::ostream & operator<<(std::ostream & fileToWrite, simol::Vector<ScalarType,WrappedLibrary> const & vectorToRead)
{
  for (size_t index = 0; index < vectorToRead.size(); ++index)
    fileToWrite << vectorToRead(index) << " ";
  return fileToWrite;
}

#include "Vector.ipp"

#endif
