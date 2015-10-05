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
}

namespace simol
{

  template<class ScalarType>
  class Vector<ScalarType,stl> : public stl<ScalarType>
  {};
  

  template<class ScalarType>
  class Vector<ScalarType,eigen>
  {
    public:
      Vector(size_t const size);
      Vector(size_t const size, ScalarType const& lambda);
      Vector(Vector<ScalarType,eigen> const& u);
      size_t size() const;
      ScalarType & operator()(size_t const index);
      ScalarType const & operator()(size_t const index) const;
      ScalarType norm() const;
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
  //std::cout << "size = " << vectorToRead.size() << std::endl;
  /*for (size_t index = 0; index < vectorToRead.size(); ++index)
  {
    std::cout << index << " / " << vectorToRead.size() << std::endl;
    fileToWrite << vectorToRead(index) << " ";
  }*/
  
  fileToWrite << vectorToRead(0) << " ";

  return fileToWrite;
}

#include "Vector.ipp"

#endif
