#ifndef VECTOR_IMPL_HPP
#define VECTOR_IMPL_HPP

namespace simol
{
  
  template<class ScalarType>
  Vector<ScalarType>::Vector(size_t size)
  :self_(size)
  {}
  
  template<class ScalarType>
  Vector<ScalarType>::Vector(size_t size, ScalarType scalar)
  :self_(size,scalar)
  {}

  template<class ScalarType>
  size_t Vector<ScalarType>::size() const
  {
    return self_.size();
  }

  template<class ScalarType>
  double & Vector<ScalarType>::operator()(size_t const index)
  {
    return self_[index];
  }

  template<class ScalarType> 
  void Vector<ScalarType>::operator*=(double scalar)
  { 
    cblas_dscal(self_.size(), scalar, &self_[0], 1);
  }

}

#endif
