#ifndef SIMOL_VECTOR_HPP
#define SIMOL_VECTOR_HPP

#include <ostream>
#include <vector>

extern "C"
{
#include <cblas.h>
}

#include "vector_forward.hpp"

namespace simol
{
  template<class ScalarType>
  std::ostream & operator<<(std::ostream & output, Vector<ScalarType> const & vector);
}

namespace simol
{

  //! \brief Vectors and vector operations
  template<class ScalarType>
  class Vector
  {
    friend std::ostream & operator<< <>(std::ostream & output, Vector<ScalarType> const & vector);
    
    public:
      //! \brief Construct a vector of a given size
      Vector(size_t size)
      :self_(size)
      {}
      //! \brief Construct a vector of a given size and assign it to a given scalar
      Vector(size_t size, ScalarType scalar)
      :self_(size,scalar)
      {}
    public:
      //! \brief Access size
      size_t size() const
      { return self_.size(); }
    public:
      double & operator()(size_t const index)
      { return self_[index]; }
    public:
      //! \brief Scaling operation
      void operator*=(double scalar)
      { cblas_dscal(self_.size(), scalar, &self_[0], 1); }
    private:
      std::vector<ScalarType> self_;

  };

  template<class ScalarType>
  std::ostream & operator<<(std::ostream & output, Vector<ScalarType> const & vector)
  {
    for (size_t index = 0; index < vector.size(); ++index)
      output << vector.self_[index] << " ";

    return output;
  }

}



#endif
