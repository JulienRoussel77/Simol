#ifndef VECTORWRAPPER_HPP
#define VECTORWRAPPER_HPP

#include "stl.hpp"

#include "Vector_fwd.hpp"

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
      size_t const & size() const;
      ScalarType & operator()(size_t const index);
      ScalarType const & operator()(size_t const index) const;

    private:
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

#include "Vector_impl.hpp"

#endif
