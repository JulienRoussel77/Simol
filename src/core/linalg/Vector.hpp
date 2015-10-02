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
      size_t size() const;
      ScalarType & operator()(size_t const index);
      ScalarType const & operator()(size_t const index) const;
      ScalarType norm() const;
      Vector<ScalarType,eigen>& operator+=(Vector<ScalarType,eigen> const& v);
      Vector<ScalarType,eigen>& operator-=(Vector<ScalarType,eigen> const& v);
      Vector<ScalarType,eigen> operator*(ScalarType const& lambda) const;
      Vector<ScalarType,eigen> operator/(ScalarType const& lambda) const;
      Vector<ScalarType,eigen> operator-() const;
      Vector<ScalarType,eigen> operator+(Vector<ScalarType,eigen> const& v) const;
      Vector<ScalarType,eigen> operator-(Vector<ScalarType,eigen> const& v) const;
      //Vector<ScalarType,eigen>& operator=(ScalarType const& lambda);	

    private:
      typename eigen<ScalarType>::VectorType wrapped_;
      
  };
  

  template<class ScalarType, template<class> class WrappedLibrary>
  std::ofstream & operator<<(std::ofstream & fileToWrite, Vector<ScalarType,WrappedLibrary> const & vectorToRead)
  {
    std::cout << "size = " << vectorToRead.size() << std::endl;
    for (size_t index = 0; index < vectorToRead.size(); ++index)
      fileToWrite << vectorToRead(index) << " ";
  
    return fileToWrite;
  }
  

  Vector<double,eigen> operator*(double const& lambda, Vector<double,eigen> const& v);
  
  


}

#include "Vector.ipp"

#endif
