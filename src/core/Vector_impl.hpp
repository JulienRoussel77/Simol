#ifndef VECTOR_IMPL_HPP
#define VECTOR_IMPL_HPP

namespace simol
{
  template<class ScalarType> inline
  Vector<ScalarType,eigen>::Vector(size_t const size)
  :wrapped_(size)
  {}
      
  template<class ScalarType> inline
  ScalarType & Vector<ScalarType,eigen>::operator()(size_t const index)
  { return wrapped_(index); }
  
  template<class ScalarType> inline
  ScalarType const & Vector<ScalarType,eigen>::operator()(size_t const index) const
  { return wrapped_(index); }
}

#endif
