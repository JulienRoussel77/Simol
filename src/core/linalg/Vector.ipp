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
  
  template<class ScalarType> inline
  Vector<ScalarType,eigen>::Vector(size_t const size, ScalarType const& lambda):wrapped_(size)
  {
   for (size_t i = 0; i<size; i++)
     wrapped_(i) = lambda;
  }
  
  template<class ScalarType> inline
  Vector<ScalarType,eigen>::Vector(Vector<ScalarType,eigen> const& u)
  :wrapped_(u.wrapped_)
  {
  }
  
  //=====================
  // ACCESSORS / MUTATORS
  //=====================

  template<class ScalarType> inline
  size_t Vector<ScalarType,eigen>::size() const
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
  
  template<class ScalarType> inline
  Vector<ScalarType,eigen>& Vector<ScalarType,eigen>::fill(ScalarType const& lambda)
  { 
    wrapped_.fill(lambda);
    return *this;
  }
  
  template<class ScalarType> inline
  ScalarType Vector<ScalarType,eigen>::dot(Vector<ScalarType,eigen> const& v) const
  {
   return wrapped_.dot(v.wrapped_); 
  }
  
  /*double dot(Vector<double,eigen> const& u, Vector<double,eigen> const& v)
  {
    return u.dot(v);
  }*/
  
  //======================
  // Operators
  //======================
  
  
  template<class ScalarType> inline
  Vector<ScalarType,eigen>& Vector<ScalarType,eigen>::operator+=(Vector<ScalarType,eigen> const& u)
  { 
    wrapped_ += u.wrapped_;
    return *this;
  }
  
  template<class ScalarType> inline
  Vector<ScalarType,eigen>& Vector<ScalarType,eigen>::operator-=(Vector<ScalarType,eigen> const& u)
  { 
    wrapped_ -= u.wrapped_;
    return *this;
  }
  
  template<class ScalarType> inline
  Vector<ScalarType,eigen>& Vector<ScalarType,eigen>::operator*=(ScalarType const& lambda)
  { 
    wrapped_ *= lambda;
    return *this;
  }
  
  template<class ScalarType> inline
  Vector<ScalarType,eigen>& Vector<ScalarType,eigen>::operator/=(ScalarType const& lambda)
  { 
    wrapped_ /= lambda;
    return *this;
  }
  
  template<class ScalarType> inline
  Vector<ScalarType,eigen> Vector<ScalarType,eigen>::operator*(ScalarType const& lambda) const
  {
    Vector<ScalarType> u(*this);
    u.wrapped_ *= lambda;
    return u;
  }
  
  template<class ScalarType> inline
  Vector<ScalarType,eigen> Vector<ScalarType,eigen>::operator/(ScalarType const& lambda) const
  {
    Vector<ScalarType> u(*this);
    u.wrapped_ /= lambda;
    return u;
  }
  
  template<class ScalarType> inline
  Vector<ScalarType,eigen> Vector<ScalarType,eigen>::operator-() const
  {
    Vector<ScalarType> u(*this);
    u.wrapped_ *= -1;
    return u;
  }
  
  template<class ScalarType> inline
  Vector<ScalarType,eigen> Vector<ScalarType,eigen>::operator+(Vector<ScalarType,eigen> const& u) const
  {
    //return wrapped_ + v.wrapped_;
    //return *this;
    return Vector<ScalarType,eigen>(*this) += u;
  }
  
  template<class ScalarType> inline
  Vector<ScalarType,eigen> Vector<ScalarType,eigen>::operator-(Vector<ScalarType,eigen> const& u) const
  {
    //return wrapped_ - v.wrapped_;
    return Vector<ScalarType,eigen>(*this) -= u;
  }
  
  
  
  // Free functions
  

  
  


}

#endif
