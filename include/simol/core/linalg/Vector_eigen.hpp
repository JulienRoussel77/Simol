#ifndef SIMOL_VECTOR_EIGEN_HPP
#define SIMOL_VECTOR_EIGEN_HPP

#include <iterator>
#include <fstream>

namespace simol
{

  template<typename Scalar>
  std::ostream & operator<<(std::ostream & fileToWrite, simol::Vector<Scalar, eigen> const & vectorToRead);

  template<class Scalar>
  class Vector<Scalar, eigen>
  {
      friend std::ostream & operator<< <>(std::ostream & fileToWrite, simol::Vector<Scalar, eigen> const & vectorToRead);

    public:

      static Vector Zero(std::size_t length);

      explicit Vector(size_t const size = 0);
      explicit Vector(size_t const size, Scalar const& lambda);
      explicit Vector(std::string const & filename);

      Vector(Vector const & vector) = default;
      Vector(typename eigen::Vector<Scalar> const & wrappedVector);

      std::size_t size() const;
      std::size_t min_index() const;

      std::vector<size_t> indices_of_smallest(size_t const number_of_indices);

      Scalar min() const;
      Scalar max() const;
      Scalar norm() const;
      Scalar & operator()(size_t const index);
      Scalar const & operator()(size_t const index) const;

      Vector sort() const;
      Vector subvector(std::size_t start, std::size_t length) const;
      Vector & fill(Scalar const& lambda);


      Vector<Scalar, eigen>& operator+=(Vector<Scalar, eigen> const& v);
      Vector<Scalar, eigen>& operator-=(Vector<Scalar, eigen> const& v);
      Vector<Scalar, eigen>& operator*=(Scalar const& lambda);
      Vector<Scalar, eigen>& operator/=(Scalar const& lambda);
      Vector<Scalar, eigen> operator*(Scalar const& lambda) const;
      Vector<Scalar, eigen> operator/(Scalar const& lambda) const;
      Vector<Scalar, eigen> operator-() const;
      Vector<Scalar, eigen> operator+(Vector<Scalar, eigen> const& v) const;
      Vector<Scalar, eigen> operator-(Vector<Scalar, eigen> const& v) const;


    public:
      typedef typename eigen::Vector<Scalar> WrappedType;
      typename eigen::Vector<Scalar> wrapped_;

  };

  Vector<double, eigen> operator*(double const& lambda, Vector<double, eigen> const& v);
  Vector<double, eigen> piecewiseDivision(Vector<double, eigen> const& u, Vector<double, eigen> const& v);
  Vector<double, eigen> piecewiseDivision(Vector<double, eigen> const& u, Vector<size_t, eigen> const& v);

  //! Returns a null vector
  template<typename Scalar>
  Vector<Scalar, eigen> Vector<Scalar, eigen>::Zero(std::size_t length)
  { return eigen::Zero<Scalar>(length); }

  //! Returns the size
  template<typename Scalar> inline
  std::size_t Vector<Scalar, eigen>::size() const
  { return eigen::size(wrapped_); }

  //! Returns the maximum coefficient
  template<typename Scalar> inline
  Scalar Vector<Scalar, eigen>::max() const
  { return eigen::max(wrapped_); }

  //! Returns the minimum coefficient
  template<typename Scalar> inline
  Scalar Vector<Scalar, eigen>::min() const
  { return eigen::min(wrapped_); }

  //! Returns the index of the minimum coefficient
  template<typename Scalar> inline
  std::size_t Vector<Scalar, eigen>::min_index() const
  { return eigen::min_index(wrapped_); };

  //! Returns a sorted copy
  template<typename Scalar> inline
  Vector<Scalar, eigen> Vector<Scalar, eigen>::sort() const
  { return eigen::sort(wrapped_); }

  //! Returns a subvector
  template<class Scalar> inline
  Vector<Scalar, eigen> Vector<Scalar, eigen>::subvector(std::size_t start, std::size_t length) const
  { return Vector<Scalar, eigen>(eigen::subvector(wrapped_, start, length)); }

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

  template<class Scalar> inline
  Vector<Scalar, eigen>::Vector(typename eigen::Vector<Scalar> const & wrappedVector)
    : wrapped_(wrappedVector)
  {}

  template<class Scalar> inline
  Vector<Scalar, eigen>::Vector(size_t const size)
    : wrapped_(size)
  {}

  template<class Scalar> inline
  Vector<Scalar, eigen>::Vector(size_t const size, Scalar const& lambda)
    : wrapped_(size)
  {
    for (size_t i = 0; i < size; i++)
      wrapped_(i) = lambda;
  }

  template<class Scalar> inline
  Vector<Scalar, eigen>::Vector(std::string const & filename):
    wrapped_(0)
  {
    std::ifstream flow(filename);
    //Reads the number of lines
    std::string line;
    int nbOfLines = 0;
    for (nbOfLines = 0; std::getline(flow, line); ++nbOfLines) {}
    wrapped_ = eigen::Vector<Scalar>(nbOfLines);
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
  template<class Scalar>
  std::vector<size_t> Vector<Scalar, eigen>::indices_of_smallest(size_t const number_of_indices)
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
  { return eigen::norm(wrapped_); }

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
  template<class Scalar> inline
  Vector<Scalar, eigen>& Vector<Scalar, eigen>::operator+=(Vector<Scalar, eigen> const& u)
  {
    wrapped_ += u.wrapped_;
    return *this;
  }

  //! substraction of another vector
  template<class Scalar> inline
  Vector<Scalar, eigen>& Vector<Scalar, eigen>::operator-=(Vector<Scalar, eigen> const& u)
  {
    wrapped_ -= u.wrapped_;
    return *this;
  }

  //! multiplication by a scalar
  template<class Scalar> inline
  Vector<Scalar, eigen>& Vector<Scalar, eigen>::operator*=(Scalar const& lambda)
  {
    wrapped_ *= lambda;
    return *this;
  }

  //! division by a scalar
  template<class Scalar> inline
  Vector<Scalar, eigen>& Vector<Scalar, eigen>::operator/=(Scalar const& lambda)
  {
    wrapped_ /= lambda;
    return *this;
  }

  template<class Scalar> inline
  Vector<Scalar, eigen> Vector<Scalar, eigen>::operator*(Scalar const& lambda) const
  {
    Vector<Scalar> u(*this);
    u.wrapped_ *= lambda;
    return u;
  }

  template<class Scalar> inline
  Vector<Scalar, eigen> Vector<Scalar, eigen>::operator/(Scalar const& lambda) const
  {
    Vector<Scalar> u(*this);
    u.wrapped_ /= lambda;
    return u;
  }

  template<class Scalar> inline
  Vector<Scalar, eigen> Vector<Scalar, eigen>::operator-() const
  {
    Vector<Scalar> u(*this);
    u.wrapped_ *= -1;
    return u;
  }

  template<class Scalar> inline
  Vector<Scalar, eigen> Vector<Scalar, eigen>::operator+(Vector<Scalar, eigen> const& u) const
  { return Vector<Scalar, eigen>(*this) += u; }

  template<class Scalar> inline
  Vector<Scalar, eigen> Vector<Scalar, eigen>::operator-(Vector<Scalar, eigen> const& u) const
  { return Vector<Scalar, eigen>(*this) -= u; }


  template<class Scalar>
  Vector<double, eigen> piecewiseDivision(Vector<double, eigen> const& u, Vector<Scalar, eigen> const& v)
  {
    if (u.size() != v.size())
      throw std::invalid_argument("Can only divide vectors of same size !");
    Vector<double, eigen> w(u.size());
    for (size_t i = 0; i < u.size(); i++)
      w(i) = u(i) / v(i);
    return w;
  }

  //! Print a vector into a file
  template<class Scalar>
  std::ostream & operator<<(std::ostream & fileToWrite, simol::Vector<Scalar, eigen> const & vectorToRead)
  {
    for (size_t index = 0; index < vectorToRead.size(); ++index)
      fileToWrite << vectorToRead(index) << " ";
    return fileToWrite;
  }


}

#endif
