#ifndef VECTOR_IMPL_HPP
#define VECTOR_IMPL_HPP

namespace simol
{
  
  template<class ScalarType, template<class> class WrappingPolicy>
  Vector<ScalarType,WrappingPolicy>::Vector(size_t size)
  :self_(size)
  {}
  
  template<class ScalarType, template<class> class WrappingPolicy>
  Vector<ScalarType,WrappingPolicy>::Vector(size_t size, ScalarType scalar)
  :self_(size,scalar)
  {}

  template<class ScalarType, template<class> class WrappingPolicy>
  size_t Vector<ScalarType,WrappingPolicy>::size() const
  {
    return self_.size();
  }

  template<class ScalarType, template<class> class WrappingPolicy>
  double & Vector<ScalarType,WrappingPolicy>::operator()(size_t const index)
  {
    return self_[index];
  }

  template<class ScalarType, template<class> class WrappingPolicy> 
  void Vector<ScalarType,WrappingPolicy>::operator*=(double scalar)
  { 
    cblas_dscal(self_.size(), scalar, &self_[0], 1);
  }

}

#endif
