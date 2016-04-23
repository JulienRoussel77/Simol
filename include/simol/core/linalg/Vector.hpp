#ifndef SIMOL_VECTOR_HPP
#define SIMOL_VECTOR_HPP

#include "eigen.hpp"

namespace simol
{
  template<typename Scalar, typename Wrapper>
  class Vector;

  template<class Scalar, typename Wrapper>
  std::ostream & operator<<(std::ostream & output, 
                            Vector<Scalar, Wrapper> const & vector);
  
  template<typename Scalar, typename Wrapper>
  Scalar inner_product(Vector<Scalar, Wrapper> const & lhs, 
                       Vector<Scalar, Wrapper> const & rhs);
  
  template<typename Scalar, typename Wrapper>
  Vector<Scalar, Wrapper> operator*(Scalar const lhs, 
                                    Vector<Scalar, Wrapper> const & rhs);
  
  template<typename Scalar, typename Wrapper>
  Vector<Scalar, Wrapper> operator-(Vector<Scalar, Wrapper> const & lhs, 
                                    Vector<Scalar, Wrapper> const & rhs);
  
  template<typename Scalar = double, typename Wrapper = eigen>
  class Vector
  {
    
    friend std::ostream & operator<<<>(std::ostream & output, Vector const & vector);
    friend Scalar inner_product<>(Vector const & lhs, Vector const & rhs);
    friend Vector operator*<>(Scalar const scalar, Vector const & rhs);
    friend Vector operator-<>(Vector const & lhs, Vector const & rhs);
    
    public:

      static Vector Zero(std::size_t const length);
      
      Vector(std::size_t const length = 0);
      Vector(std::size_t const length, Scalar const scalar);
      Vector(typename eigen::Vector<Scalar> const & vector);
      
      void fill(Scalar const scalar);


      std::size_t size() const;
      std::size_t min_index() const;

      std::vector<size_t> indices_of_smallest(size_t const number_of_indices);
      
      Scalar min() const;
      Scalar max() const;
      Scalar norm() const;
      Scalar & operator()(size_t const index);
      Scalar const & operator()(size_t const index) const;

      Vector sort() const;
      Vector subvector(std::size_t const start, std::size_t const length) const;
      Vector operator+(Vector const & vector) const;
      Vector operator*(Scalar const scalar) const;
      Vector operator/(Scalar const scalar) const;
      Vector operator-() const;
      Vector & operator*=(Scalar const scalar);
      Vector & operator/=(Scalar const scalar);
      Vector & operator+=(Vector const & vector);
      Vector & operator-=(Vector const & vector);

    public:

      typename Wrapper::template Vector<Scalar> wrapped_;

  };

  //! Returns a null vector
  template<typename Scalar, typename Wrapper> inline
  Vector<Scalar, Wrapper> Vector<Scalar, Wrapper>::Zero(std::size_t const length)
  { return Vector<Scalar, Wrapper>( Wrapper::template Zero<Scalar>(length) ); }

  //! Fill a vector with some scalar
  template<typename Scalar, typename Wrapper> inline
  void Vector<Scalar, Wrapper>::fill(Scalar const scalar)
  { Wrapper::fill(wrapped_, scalar); }

  //! Construction of a vector from a length
  template<typename Scalar, typename Wrapper> inline
  Vector<Scalar, Wrapper>::Vector(std::size_t const length)
  : wrapped_(length)
  {}

  //! Construction of a vector from a length
  template<typename Scalar, typename Wrapper> inline
  Vector<Scalar, Wrapper>::Vector(std::size_t const length, Scalar const scalar)
  : wrapped_(length)
  { Wrapper::fill(wrapped_, scalar); }

  //! Construction of a vector from a wrapped vector
  template<typename Scalar, typename Wrapper> inline
  Vector<Scalar, Wrapper>::Vector(typename eigen::Vector<Scalar> const & vector)
  : wrapped_(vector)
  {}

  // Returns the size
  template<typename Scalar, typename Wrapper> inline
  std::size_t Vector<Scalar, Wrapper>::size() const
  { return Wrapper::size(wrapped_); }

  //! Returns the index of the minimum coefficient
  template<typename Scalar, typename Wrapper> inline
  std::size_t Vector<Scalar, Wrapper>::min_index() const
  { return Wrapper::min_index(wrapped_); };

  //! Returns the indices of the first smallest coefficients
  template<typename Scalar, typename Wrapper> inline
  std::vector<size_t> Vector<Scalar, Wrapper>::indices_of_smallest(size_t const number_of_indices)
  {
    std::vector<size_t> indices(size());
    for(size_t index = 0; index < size(); ++index)
      indices[index] = index;

    std::sort(indices.begin(),
              indices.end(),
              [&](size_t const & left_index,
                  size_t const & right_index)
    { return (wrapped_[left_index] < wrapped_[right_index]); }
             );

    return std::vector<size_t>(indices.begin(), indices.begin() + number_of_indices);

  }

  //! Returns the maximum coefficient
  template<typename Scalar, typename Wrapper> inline
  Scalar Vector<Scalar, Wrapper>::max() const
  { return Wrapper::max(wrapped_); }

  //! Returns the minimum coefficient
  template<typename Scalar, typename Wrapper> inline
  Scalar Vector<Scalar, Wrapper>::min() const
  { return Wrapper::min(wrapped_); }

  //! Returns the Euclidean norm
  template<typename Scalar, typename Wrapper> inline
  Scalar Vector<Scalar, Wrapper>::norm() const
  { return Wrapper::norm(wrapped_); }

  // Returns a mutable coefficient
  template<typename Scalar, typename Wrapper> inline
  Scalar & Vector<Scalar, Wrapper>::operator()(size_t const index)
  { return wrapped_(index); }

  // Returns an immutable coefficient
  template<typename Scalar, typename Wrapper> inline
  Scalar const & Vector<Scalar, Wrapper>::operator()(size_t const index) const
  { return wrapped_(index); }

  //! Returns a sorted copy
  template<typename Scalar, typename Wrapper> inline
  Vector<Scalar, Wrapper> Vector<Scalar, Wrapper>::sort() const
  { return Wrapper::sort(wrapped_); }

  //! Returns a subvector
  template<typename Scalar, typename Wrapper> inline
  Vector<Scalar, Wrapper> Vector<Scalar, Wrapper>::subvector(std::size_t const start, std::size_t const length) const
  { return Vector<Scalar, Wrapper>(Wrapper::subvector(wrapped_, start, length)); }

  // Sum of vectors
  template<typename Scalar, typename Wrapper> inline
  Vector<Scalar, Wrapper> Vector<Scalar, Wrapper>::operator+(Vector<Scalar, Wrapper> const& u) const
  { return Vector<Scalar, Wrapper>(*this) += u; }

  //! Multiplication by a scalar
  template<typename Scalar, typename Wrapper> inline
  Vector<Scalar, Wrapper> Vector<Scalar, Wrapper>::operator*(Scalar const scalar) const
  {
    Vector<Scalar, Wrapper> u(*this);
    u.wrapped_ *= scalar;
    return u;
  }

  //! Division by a scalar
  template<typename Scalar, typename Wrapper> inline
  Vector<Scalar, Wrapper> Vector<Scalar, Wrapper>::operator/(Scalar const scalar) const
  {
    Vector<Scalar, Wrapper> u(*this);
    u.wrapped_ /= scalar;
    return u;
  }

  //! Change the sign
  template<typename Scalar, typename Wrapper> inline
  Vector<Scalar, Wrapper> Vector<Scalar, Wrapper>::operator-() const
  {
    Vector<Scalar, Wrapper> u(*this);
    u.wrapped_ *= -1;
    return u;
  }

  //! Multiply by a scalar
  template<typename Scalar, typename Wrapper> inline
  Vector<Scalar, Wrapper> & Vector<Scalar, Wrapper>::operator*=(Scalar const scalar)
  { 
    Wrapper::multiply(wrapped_, scalar);
    return *this;
  }

  //! Divide by a scalar
  template<typename Scalar, typename Wrapper> inline
  Vector<Scalar, Wrapper> & Vector<Scalar, Wrapper>::operator/=(Scalar const scalar)
  { 
    wrapped_ /= scalar;
    return *this;
  }

  //! Add a vector
  template<typename Scalar, typename Wrapper> inline
  Vector<Scalar, Wrapper> & Vector<Scalar, Wrapper>::operator+=(Vector<Scalar, Wrapper> const & vector)
  { 
    Wrapper::add(wrapped_, vector.wrapped_);
    return *this;
  }

  //! Substract a vector
  template<typename Scalar, typename Wrapper> inline
  Vector<Scalar, Wrapper> & Vector<Scalar, Wrapper>::operator-=(Vector<Scalar, Wrapper> const & vector)
  { 
    wrapped_ -= vector.wrapped_;
    return *this;
  }

  //! Print a vector into a file
  template<class Scalar, typename Wrapper> inline
  std::ostream & operator<<(std::ostream & output, 
                            Vector<Scalar, Wrapper> const & vector)
  {
    for (size_t index = 0; index < vector.size(); ++index)
      output << vector(index) << " ";
    return output;
  }

  //! Vector inner product
  template<typename Scalar, typename Wrapper> inline
  Scalar inner_product(Vector<Scalar, Wrapper> const & lhs, 
                       Vector<Scalar, Wrapper> const & rhs)
  { return Wrapper::inner_product(lhs.wrapped_, rhs.wrapped_); }

  //! Product between a scalar and a vector
  template<typename Scalar, typename Wrapper> inline
  Vector<Scalar, Wrapper> operator*(Scalar const lhs, 
                                    Vector<Scalar, Wrapper> const & rhs)
  { return Vector<Scalar, Wrapper>(lhs * rhs); }

  //! Vector difference
  template<typename Scalar, typename Wrapper> inline
  Vector<Scalar, Wrapper> operator-(Vector<Scalar, Wrapper> const & lhs, 
                                    Vector<Scalar, Wrapper> const & rhs)
  { return Vector<Scalar, Wrapper>(lhs - rhs); }

  template<typename Scalar, typename Wrapper> inline
  Vector<double, eigen> piecewiseDivision(Vector<Scalar, Wrapper> const & lhs, 
                                          Vector<int, Wrapper> const & rhs)
  {
    if (lhs.size() != rhs.size())
      throw std::invalid_argument("Can only divide vectors of same size !");
    Vector<double, Wrapper> division(lhs.size());
    for (size_t i = 0; i < lhs.size(); i++)
      division(i) = lhs(i) / rhs(i);
    return division;
  }


}


#endif
