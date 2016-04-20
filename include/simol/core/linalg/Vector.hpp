#ifndef SIMOL_VECTOR_HPP
#define SIMOL_VECTOR_HPP

#include "VectorInterface.hpp"

#include "eigen.hpp"

namespace simol
{
  template<typename Scalar = double, template<class> class Wrapper = eigen>
  class Vector
  {
    public:

      std::size_t size() const;
      std::size_t min_index() const;

      Scalar min() const;
      Scalar max() const;
      Scalar norm() const;

      Vector sort() const;

    private:

      Wrapper<Scalar> wrapped_;

  };

  template<typename Scalar, template<class> class Wrapper> inline
  std::size_t Vector<Scalar, Wrapper>::size() const
  { return Wrapper<Scalar>::size(wrapped_); }


}

#include "Vector_eigen.hpp"

#endif
