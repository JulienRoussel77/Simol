#ifndef VECTORWRAPPER_HPP
#define VECTORWRAPPER_HPP

#include "stl.hpp"
#include "eigen.hpp"
#include <iostream>

namespace simol
{
  template<class ScalarType=double, template<class> class WrappedLibrary=eigen>
  class Vector;

  typedef Vector<double> dvec;

}

namespace simol
{

  template<class ScalarType>
  class Vector<ScalarType,stl> : public stl<ScalarType>
  {};


  template<class ScalarType>
  class Vector<ScalarType,eigen>
  {
    //friend double dot(Vector<double,eigen> const& u, Vector<double,eigen> const& v);
    public:
      Vector(size_t const size=0);
      Vector(size_t const size, ScalarType const& lambda);
      Vector(Vector<ScalarType,eigen> const& u);
      Vector(typename eigen<ScalarType>::VectorType const & wrappedVector);
      size_t size() const;
      ScalarType & operator()(size_t const index);
      ScalarType const & operator()(size_t const index) const;
      ScalarType norm() const;
      Vector<ScalarType,eigen>& fill(ScalarType const& lambda);
      ScalarType dot(Vector<ScalarType,eigen> const& v) const;
      ScalarType min() const;
      size_t min_index() const;
      Vector<ScalarType, eigen> sort() const;
      std::vector<size_t> smallest_indices(size_t const number_of_indices);

      ScalarType max() const
      { return wrapped_.maxCoeff(); }
      Vector<ScalarType,eigen>& operator+=(Vector<ScalarType,eigen> const& v);
      Vector<ScalarType,eigen>& operator-=(Vector<ScalarType,eigen> const& v);
      Vector<ScalarType,eigen>& operator*=(ScalarType const& lambda);
      Vector<ScalarType,eigen>& operator/=(ScalarType const& lambda);
      Vector<ScalarType,eigen> operator*(ScalarType const& lambda) const;
      Vector<ScalarType,eigen> operator/(ScalarType const& lambda) const;
      Vector<ScalarType,eigen> operator-() const;
      Vector<ScalarType,eigen> operator+(Vector<ScalarType,eigen> const& v) const;
      Vector<ScalarType,eigen> operator-(Vector<ScalarType,eigen> const& v) const;


    public:
      typename eigen<ScalarType>::VectorType wrapped_;

  };

  /*template<class ScalarType, template<class> class WrappedLibrary>
  std::ofstream & operator<<(std::ofstream & fileToWrite, Vector<ScalarType,WrappedLibrary> const & vectorToRead)
  {
    for (size_t index = 0; index < vectorToRead.size(); ++index)
      fileToWrite << vectorToRead(index) << " ";

    return fileToWrite;
  }*/

  Vector<double,eigen> operator*(double const& lambda, Vector<double,eigen> const& v);
  double dot(Vector<double,eigen> const& u, Vector<double,eigen> const& v);

  //=============
  // CONSTRUCTORS
  //=============

  template<class ScalarType>
  inline
  Vector<ScalarType, eigen>::Vector(typename eigen<ScalarType>::VectorType const & wrappedVector)
  :wrapped_(wrappedVector)
  {}

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
  {}

  //=====================
  // ACCESSORS / MUTATORS
  //=====================

  template<class ScalarType>
  inline
  size_t
  Vector<ScalarType,eigen>::size() const
  { return wrapped_.size(); }

  template<class ScalarType>
  inline
  ScalarType &
  Vector<ScalarType,eigen>::operator()(size_t const index)
  { return wrapped_(index); }

  template<class ScalarType>
  inline
  ScalarType const &
  Vector<ScalarType,eigen>::operator()(size_t const index) const
  { return wrapped_(index); }

  //======================
  // Utils
  //======================

  template<class ScalarType>
  inline
  ScalarType
  Vector<ScalarType, eigen>::min() const
  { return wrapped_.minCoeff(); }

  template<class ScalarType>
  inline
  Vector<ScalarType, eigen>
  Vector<ScalarType, eigen>::sort() const
  {
      Vector<ScalarType, eigen> to_be_sorted = *this;
      std::sort( to_be_sorted.wrapped_.data(), to_be_sorted.wrapped_.data() + size() );
      return to_be_sorted;
  }

  template<class ScalarType>
  inline
  size_t
  Vector<ScalarType, eigen>::min_index() const
  {
      size_t index;
      wrapped_.minCoeff(&index);
      return index;
  }

  template<class ScalarType>
  inline
  std::vector<size_t>
  Vector<ScalarType, eigen>::smallest_indices(size_t const number_of_indices)
  {
    std::vector<size_t> indices(size());
    for(size_t index=0; index<size(); ++index)
        indices[index] = index;

    std::sort(indices.begin(),
              indices.end(),
              [&](size_t const & left_index,
                  size_t const & right_index)
              { return (wrapped_[left_index] < wrapped_[right_index]); }
             );

    return std::vector<size_t>(indices.begin(), indices.begin() + number_of_indices);

  }


  template<class ScalarType>
  inline
  ScalarType
  Vector<ScalarType,eigen>::norm() const
  { return wrapped_.norm(); }

  template<class ScalarType>
  inline
  Vector<ScalarType,eigen> &
  Vector<ScalarType,eigen>::fill(ScalarType const& lambda)
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

}

template<class ScalarType, template<class> class WrappedLibrary>
std::ostream & operator<<(std::ostream & fileToWrite, simol::Vector<ScalarType,WrappedLibrary> const & vectorToRead)
{
  for (size_t index = 0; index < vectorToRead.size(); ++index)
    fileToWrite << vectorToRead(index) << " ";
  return fileToWrite;
}


#endif
