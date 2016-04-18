#ifndef SIMOL_VECTOR_INTERFACE_HPP
#define SIMOL_VECTOR_INTERFACE_HPP

#include <cstddef>

namespace simol
{
  template<typename Scalar, template<class> class Library>
  class Vector;

  template<typename Scalar, template<class> class Library>
  class VectorInterface
  {
    public:
      Scalar & operator()(std::size_t const index);

    protected:
      ~VectorInterface() = default;
  };

  template<typename Scalar, template<class> class Library>
  Scalar & VectorInterface<Scalar, Library>::operator()(std::size_t const index)
  { return static_cast< Vector<Scalar, Library>* >(this)->operator()(index); }


}


#endif
