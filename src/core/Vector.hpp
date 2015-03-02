#ifndef SIMOL_VECTOR_HPP
#define SIMOL_VECTOR_HPP

#include <ostream>
#include <vector>

extern "C"
{
#include <cblas.h>
}

#include "Vector_fwd.hpp"

namespace simol
{
  template<class ScalarType, template<class> class WrappingPolicy>
  std::ostream & operator<<(std::ostream & output, Vector<ScalarType,WrappingPolicy> const & vector);
}

namespace simol
{

  template<class ScalarType> class Dummy {};

  //! \brief Vectors and vector operations
  template<class ScalarType, template<class> class WrappingPolicy>
  class Vector
  {
    friend std::ostream & operator<< <>(std::ostream & output, Vector<ScalarType,WrappingPolicy> const & vector);
    
    public:
      //! \brief Construct a vector of a given size
      Vector(size_t size);
      
      //! \brief Construct a vector of a given size and assign it to a given scalar
      Vector(size_t size, ScalarType scalar);
    public:
      //! \brief Access size
      size_t size() const;
    public:
      double & operator()(size_t const index);
    public:
      //! \brief Scaling operation
      void operator*=(double scalar);
    private:
      std::vector<ScalarType> self_;

  };

}

#include "Vector_impl.hpp"

namespace simol
{
  template<class ScalarType, template<class> class WrappingPolicy>
  std::ostream & operator<<(std::ostream & output, Vector<ScalarType,WrappingPolicy> const & vector)
  {
    for (size_t index = 0; index < vector.size(); ++index)
      output << vector.self_[index] << " ";

    return output;
  }

}


#endif
