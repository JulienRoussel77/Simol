#ifndef SIMOL_VECTOR_EIGEN_HPP
#define SIMOL_VECTOR_EIGEN_HPP

#include <iterator>
#include <fstream>

namespace simol
{

  template<typename ScalarType>
  std::ostream & operator<<(std::ostream & fileToWrite, simol::Vector<ScalarType, eigen> const & vectorToRead);

  template<class ScalarType>
  class Vector<ScalarType, eigen> : public VectorInterface<ScalarType, eigen>
  {
      friend std::ostream & operator<< <>(std::ostream & fileToWrite, simol::Vector<ScalarType, eigen> const & vectorToRead);

    public:

      static Vector Zero(std::size_t length);
      
      explicit Vector(size_t const size = 0);
      explicit Vector(size_t const size, ScalarType const& lambda);
      explicit Vector(std::string const & filename);

      Vector(Vector<ScalarType, eigen> const& u);
      Vector(typename eigen<ScalarType>::Vector const & wrappedVector);

      std::size_t size() const;
      std::size_t min_index() const;
      
      std::vector<size_t> indices_of_smallest(size_t const number_of_indices);
      
      ScalarType min() const;
      ScalarType max() const;
      ScalarType norm() const;
      ScalarType & operator()(size_t const index);
      ScalarType const & operator()(size_t const index) const;
      
      Vector sort() const;
      Vector subvector(std::size_t start, std::size_t length) const;
      Vector & fill(ScalarType const& lambda);


      Vector<ScalarType, eigen>& operator+=(Vector<ScalarType, eigen> const& v);
      Vector<ScalarType, eigen>& operator-=(Vector<ScalarType, eigen> const& v);
      Vector<ScalarType, eigen>& operator*=(ScalarType const& lambda);
      Vector<ScalarType, eigen>& operator/=(ScalarType const& lambda);
      Vector<ScalarType, eigen> operator*(ScalarType const& lambda) const;
      Vector<ScalarType, eigen> operator/(ScalarType const& lambda) const;
      Vector<ScalarType, eigen> operator-() const;
      Vector<ScalarType, eigen> operator+(Vector<ScalarType, eigen> const& v) const;
      Vector<ScalarType, eigen> operator-(Vector<ScalarType, eigen> const& v) const;


    public:
      typedef typename eigen<ScalarType>::Vector WrappedType;
      typename eigen<ScalarType>::Vector wrapped_;

  };

  Vector<double, eigen> operator*(double const& lambda, Vector<double, eigen> const& v);
  Vector<double, eigen> piecewiseDivision(Vector<double, eigen> const& u, Vector<double, eigen> const& v);
  Vector<double, eigen> piecewiseDivision(Vector<double, eigen> const& u, Vector<size_t, eigen> const& v);

  //! Returns a null vector
  template<typename Scalar>
  Vector<Scalar, eigen> Vector<Scalar, eigen>::Zero(std::size_t length)
  { return eigen<Scalar>::Zero(length); }

  //! Returns the size
  template<typename Scalar> inline
  std::size_t Vector<Scalar, eigen>::size() const
  { return eigen<Scalar>::size(wrapped_); }

  //! Returns the maximum coefficient
  template<typename Scalar> inline
  Scalar Vector<Scalar, eigen>::max() const
  { return eigen<Scalar>::max(wrapped_); }

  //! Returns the minimum coefficient
  template<typename Scalar> inline
  Scalar Vector<Scalar, eigen>::min() const
  { return eigen<Scalar>::min(wrapped_); }

  //! Returns the index of the minimum coefficient
  template<typename Scalar> inline
  std::size_t Vector<Scalar, eigen>::min_index() const
  { return eigen<Scalar>::min_index(wrapped_); };

  //! Returns a sorted copy
  template<typename Scalar> inline
  Vector<Scalar, eigen> Vector<Scalar, eigen>::sort() const
  { return eigen<Scalar>::sort(wrapped_); }

  //! Returns a subvector
  template<class Scalar> inline
  Vector<Scalar, eigen> Vector<Scalar, eigen>::subvector(std::size_t start, std::size_t length) const
  { return Vector<Scalar, eigen>(eigen<Scalar>::subvector(wrapped_, start, length)); }

  //! Returns a coefficient to be modified
  template<typename Scalar> inline
  Scalar & Vector<Scalar, eigen>::operator()(size_t const index)
  { return wrapped_(index); }

  // Returns a non-modifiable coefficient 
  template<class Scalar> inline
  Scalar const & Vector<Scalar, eigen>::operator()(size_t const index) const
  { return wrapped_(index); }



  //=============
  // CONSTRUCTORS
  //=============

  template<class ScalarType> inline
  Vector<ScalarType, eigen>::Vector(typename eigen<ScalarType>::Vector const & wrappedVector)
  : wrapped_(wrappedVector)
  {}

  template<class ScalarType> inline
  Vector<ScalarType, eigen>::Vector(size_t const size)
  : wrapped_(size)
  {}

  template<class ScalarType> inline
  Vector<ScalarType, eigen>::Vector(size_t const size, ScalarType const& lambda)
  : wrapped_(size)
  {
    for (size_t i = 0; i < size; i++)
      wrapped_(i) = lambda;
  }

  template<class ScalarType> inline
  Vector<ScalarType, eigen>::Vector(Vector<ScalarType, eigen> const& u)
    : wrapped_(u.wrapped_)
  {}

  template<class ScalarType> inline
  Vector<ScalarType, eigen>::Vector(std::string const & filename):
    wrapped_(0)
  {
    std::ifstream flow(filename);
    //Reads the number of lines
    std::string line;
    int nbOfLines = 0;
    for (nbOfLines = 0; std::getline(flow, line); ++nbOfLines) {}
    wrapped_ = eigen<ScalarType>::VectorType(nbOfLines);
    //Go back to beginning of file
    flow.clear();
    flow.seekg(0, std::ios::beg);
    //Reads the coefficients
    for (size_t i = 0; i < nbOfLines; i++)
    {
      std::getline(flow, line);
      wrapped_(i) = std::stod(line);
    }
  }

  //=====================
  // ACCESSORS / MUTATORS
  //=====================

  //----- mathematical functions -----

  //! Returns the indices of the first smallest coefficients
  template<class ScalarType>
  std::vector<size_t> Vector<ScalarType, eigen>::indices_of_smallest(size_t const number_of_indices)
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

  //! Returns the Euclidean norm
  template<typename Scalar> inline
  Scalar Vector<Scalar, eigen>::norm() const
  { return eigen<Scalar>::norm(wrapped_); }

  template<typename Scalar> inline
  Vector<Scalar, eigen> & Vector<Scalar, eigen>::fill(Scalar const& lambda)
  {
    wrapped_.fill(lambda);
    return *this;
  }

  template<typename Scalar> inline
  Scalar inner_product(Vector<Scalar, eigen> const& v, Vector<Scalar, eigen> const & w)
  { return inner_product(v, w); }

  //======================
  // Operators
  //======================


  //! addition of another vector
  template<class ScalarType> inline
  Vector<ScalarType, eigen>& Vector<ScalarType, eigen>::operator+=(Vector<ScalarType, eigen> const& u)
  {
    wrapped_ += u.wrapped_;
    return *this;
  }

  //! substraction of another vector
  template<class ScalarType> inline
  Vector<ScalarType, eigen>& Vector<ScalarType, eigen>::operator-=(Vector<ScalarType, eigen> const& u)
  {
    wrapped_ -= u.wrapped_;
    return *this;
  }

  //! multiplication by a scalar
  template<class ScalarType> inline
  Vector<ScalarType, eigen>& Vector<ScalarType, eigen>::operator*=(ScalarType const& lambda)
  {
    wrapped_ *= lambda;
    return *this;
  }

  //! division by a scalar
  template<class ScalarType> inline
  Vector<ScalarType, eigen>& Vector<ScalarType, eigen>::operator/=(ScalarType const& lambda)
  {
    wrapped_ /= lambda;
    return *this;
  }

  template<class ScalarType> inline
  Vector<ScalarType, eigen> Vector<ScalarType, eigen>::operator*(ScalarType const& lambda) const
  {
    Vector<ScalarType> u(*this);
    u.wrapped_ *= lambda;
    return u;
  }

  template<class ScalarType> inline
  Vector<ScalarType, eigen> Vector<ScalarType, eigen>::operator/(ScalarType const& lambda) const
  {
    Vector<ScalarType> u(*this);
    u.wrapped_ /= lambda;
    return u;
  }

  template<class ScalarType> inline
  Vector<ScalarType, eigen> Vector<ScalarType, eigen>::operator-() const
  {
    Vector<ScalarType> u(*this);
    u.wrapped_ *= -1;
    return u;
  }

  template<class ScalarType> inline
  Vector<ScalarType, eigen> Vector<ScalarType, eigen>::operator+(Vector<ScalarType, eigen> const& u) const
  { return Vector<ScalarType, eigen>(*this) += u; }

  template<class ScalarType> inline
  Vector<ScalarType, eigen> Vector<ScalarType, eigen>::operator-(Vector<ScalarType, eigen> const& u) const
  { return Vector<ScalarType, eigen>(*this) -= u; }


  template<class ScalarType>
  Vector<double, eigen> piecewiseDivision(Vector<double, eigen> const& u, Vector<ScalarType, eigen> const& v)
  {
    if (u.size() != v.size())
      throw std::invalid_argument("Can only divide vectors of same size !");
    Vector<double, eigen> w(u.size());
    for (size_t i = 0; i < u.size(); i++)
      w(i) = u(i) / v(i);
    return w;
  }

  //! Print a vector into a file
  template<class ScalarType>
  std::ostream & operator<<(std::ostream & fileToWrite, simol::Vector<ScalarType, eigen> const & vectorToRead)
  {
    for (size_t index = 0; index < vectorToRead.size(); ++index)
      fileToWrite << vectorToRead(index) << " ";
    return fileToWrite;
  }


}

#endif
