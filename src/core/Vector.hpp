#ifndef VECTOR_HPP
#define VECTOR_HPP

#include "VectorWrapper.hpp"

namespace simol
{

  template<class ScalarType, template<class> class WrappedLibrary = eigen>
  using Vector = typename VectorWrapper<ScalarType,WrappedLibrary>::VectorType;


}
  template<class ScalarType>
  std::ofstream & operator<<(std::ofstream & fileToWrite, std::vector<ScalarType> const & vectorToRead)
  {
    for (size_t index = 0; index < vectorToRead.size(); ++index)
      fileToWrite << vectorToRead[index] << " ";
  
    return fileToWrite;
  }

#endif
