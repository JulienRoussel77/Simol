#ifndef SIMOL_VECTOR_HPP
#define SIMOL_VECTOR_HPP

#include "eigen.hpp"

namespace simol
{
  template<typename Scalar = double, typename Wrapper = eigen>
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

      typename Wrapper::template Vector<Scalar> wrapped_;

  };

  template<typename Scalar, typename Wrapper> inline
  std::size_t Vector<Scalar, Wrapper>::size() const
  { return Wrapper::size(wrapped_); }


}

#include "Vector_eigen.hpp"

#endif
