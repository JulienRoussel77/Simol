#ifndef VECTORWRAPPER_HPP
#define VECTORWRAPPER_HPP

#include "stl.hpp"
#include "eigen.hpp"

namespace simol
{

  template<class ScalarType, template<class> class WrappedLibrary = eigen>
  class Vector;
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
      Vector(typename eigen<ScalarType>::VectorType const & wrappedVector)
      :wrapped_(wrappedVector)
      {}
      size_t size() const;
      ScalarType & operator()(size_t const index);
      ScalarType const & operator()(size_t const index) const;
      
      ScalarType min() const
      { return wrapped_.minCoeff(); }
      
      ScalarType max() const
      { return wrapped_.maxCoeff(); }
      
      

    public:
      typename eigen<ScalarType>::VectorType wrapped_;
      
  };

  template<class ScalarType, template<class> class WrappedLibrary>
  std::ofstream & operator<<(std::ofstream & fileToWrite, Vector<ScalarType,WrappedLibrary> const & vectorToRead)
  {
    for (size_t index = 0; index < vectorToRead.size(); ++index)
      fileToWrite << vectorToRead(index) << " ";
  
    return fileToWrite;
  }


}

#include "Vector.ipp"

#endif
