#ifndef VECTOR_IMPL_HPP
#define VECTOR_IMPL_HPP

namespace simol
{

  //=============
  // CONSTRUCTORS
  //=============

  template<class ScalarType> inline
  Vector<ScalarType,eigen>::Vector(size_t const size)
  :wrapped_(size)
  {}
  
  //=====================
  // ACCESSORS / MUTATORS
  //=====================

  template<class ScalarType> inline
  size_t const & Vector<ScalarType,eigen>::size() const
  { return wrapped_.size(); }
  
  template<class ScalarType> inline
  ScalarType & Vector<ScalarType,eigen>::operator()(size_t const index)
  { return wrapped_(index); }
  
  template<class ScalarType> inline
  ScalarType const & Vector<ScalarType,eigen>::operator()(size_t const index) const
  { return wrapped_(index); }
  
  //======================
  // Utils
  //======================
  
  template<class ScalarType> inline
  ScalarType Vector<ScalarType,eigen>::norm() const
  { return wrapped_.norm(); }
  
  //======================
  // Operators
  //======================
  
  template<class ScalarType> inline
  Vector<ScalarType,eigen>& Vector<ScalarType,eigen>::operator+=(Vector<ScalarType,eigen> const& v)
  { 
    return *this;
  }
  
  template<class ScalarType> inline
  Vector<ScalarType,eigen>& Vector<ScalarType,eigen>::operator-=(Vector<ScalarType,eigen> const& v)
  { 
    return *this;
  }
  
  template<class ScalarType> inline
  Vector<ScalarType,eigen> Vector<ScalarType,eigen>::operator*(ScalarType const& lambda) const
  {
    //return wrapped_ * lambda;
    return *this;
  }
  
  template<class ScalarType> inline
  Vector<ScalarType,eigen> Vector<ScalarType,eigen>::operator/(ScalarType const& lambda) const
  {
    //return wrapped_ / lambda;
    return *this;
  }
  
  template<class ScalarType> inline
  Vector<ScalarType,eigen> Vector<ScalarType,eigen>::operator-() const
  {
    //return -wrapped_;
    return *this;
  }
  
  template<class ScalarType> inline
  Vector<ScalarType,eigen> Vector<ScalarType,eigen>::operator+(Vector<ScalarType,eigen> const& v) const
  {
    //return wrapped_ + v.wrapped_;
    return *this;
  }
  
  template<class ScalarType> inline
  Vector<ScalarType,eigen> Vector<ScalarType,eigen>::operator-(Vector<ScalarType,eigen> const& v) const
  {
    //return wrapped_ - v.wrapped_;
    return *this;
  }
  
  /*template<class ScalarType> inline
  Vector<ScalarType,eigen>& Vector<ScalarType,eigen>::operator=(ScalarType const& lambda)
  {
    for (int i=1; 0<size(); i++)
      wrapped_(i) = lambda;
    return *this;
  }*/
  
  
  
  // Free functions
  

  
  


}

#endif
