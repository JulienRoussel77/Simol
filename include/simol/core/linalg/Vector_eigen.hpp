#ifndef SIMOL_VECTOR_EIGEN_HPP
#define SIMOL_VECTOR_EIGEN_HPP

#include <iterator>
#include <fstream>

namespace simol
{

  template<typename ScalarType>
  std::ostream & operator<<(std::ostream & fileToWrite, simol::Vector<ScalarType,eigen> const & vectorToRead);

  template<class ScalarType>
  class Vector<ScalarType,eigen> : public VectorInterface<ScalarType, eigen>
  {
    friend std::ostream & operator<< <>(std::ostream & fileToWrite, simol::Vector<ScalarType,eigen> const & vectorToRead);
     
    public:

      explicit Vector(size_t const size=0);
      explicit Vector(size_t const size, ScalarType const& lambda);
      explicit Vector(std::string const & filename);
      
      Vector(Vector<ScalarType,eigen> const& u);
      Vector(typename eigen<ScalarType>::VectorType const & wrappedVector);
      
      size_t size() const;
      
      ScalarType & operator()(size_t const index);
      ScalarType const & operator()(size_t const index) const;
      
      ScalarType norm() const;
      Vector<ScalarType,eigen>& fill(ScalarType const& lambda);
      ScalarType dot(Vector<ScalarType,eigen> const& v) const;
      Vector<ScalarType, eigen> sort() const;
      std::vector<size_t> indices_of_smallest(size_t const number_of_indices);
      static Vector Zero(std::size_t length);

      ScalarType min() const;

      ScalarType max() const; 

      size_t index_of_minimum() const;
      
      Vector<ScalarType,eigen>& operator+=(Vector<ScalarType,eigen> const& v);
      Vector<ScalarType,eigen>& operator-=(Vector<ScalarType,eigen> const& v);
      Vector<ScalarType,eigen>& operator*=(ScalarType const& lambda);
      Vector<ScalarType,eigen>& operator/=(ScalarType const& lambda);
      Vector<ScalarType,eigen> operator*(ScalarType const& lambda) const;
      Vector<ScalarType,eigen> operator/(ScalarType const& lambda) const;
      Vector<ScalarType,eigen> operator-() const;
      Vector<ScalarType,eigen> operator+(Vector<ScalarType,eigen> const& v) const;
      Vector<ScalarType,eigen> operator-(Vector<ScalarType,eigen> const& v) const;

      Vector<ScalarType, eigen> subvec(std::size_t start, std::size_t length) const;

    public:
      typedef typename eigen<ScalarType>::VectorType WrappedType;
      typename eigen<ScalarType>::VectorType wrapped_;

  };

  Vector<double,eigen> operator*(double const& lambda, Vector<double,eigen> const& v);
  double dot(Vector<double,eigen> const& u, Vector<double,eigen> const& v);
  Vector<double,eigen> piecewiseDivision(Vector<double,eigen> const& u, Vector<double,eigen> const& v);
  Vector<double,eigen> piecewiseDivision(Vector<double,eigen> const& u, Vector<size_t,eigen> const& v);


  //----- static functions -----
  
  template<typename ScalarType>
  Vector<ScalarType, eigen> Vector<ScalarType, eigen>::Zero(std::size_t length)
  { return WrappedType(WrappedType::Zero(length)); }


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
  
  template<class ScalarType> inline
  Vector<ScalarType,eigen>::Vector(std::string const & filename):
    wrapped_(0)
  {
    std::ifstream flow(filename);
    //Reads the number of lines
    std::string line;
    int nbOfLines=0;
    for (nbOfLines = 0; std::getline(flow, line); ++nbOfLines){}
    wrapped_ = eigen<ScalarType>::VectorType(nbOfLines);
    //Go back to beginning of file
    flow.clear();
    flow.seekg(0, std::ios::beg);
    //Reads the coefficients
    for (size_t i = 0; i<nbOfLines; i++)
    {
      std::getline(flow, line);
      wrapped_(i) = std::stod(line);
    }
  }


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

  template<class ScalarType>
  Vector<ScalarType, eigen> Vector<ScalarType, eigen>::subvec(std::size_t start,
                                                              std::size_t length) const
  {
    Vector sub(length);
    for (std::size_t index = 0; index < length; ++index)
      sub.wrapped_(index) = wrapped_(index+start);
    return sub;
  }




  //----- mathematical functions -----
  
  //! Returns the maximum coefficient
  template<typename ScalarType> inline
  ScalarType Vector<ScalarType, eigen>::max() const
  { return wrapped_.maxCoeff(); }

  //! Returns the minimum coefficient
  template<class ScalarType> inline
  ScalarType Vector<ScalarType, eigen>::min() const
  { return wrapped_.minCoeff(); }
  
  //! Returns the index of the minimum coefficient
  template<class ScalarType> inline
  std::size_t Vector<ScalarType, eigen>::index_of_minimum() const
  {
      size_t index;
      wrapped_.minCoeff(&index);
      return index;
  }

  //! Returns a sorted copy
  template<class ScalarType> inline
  Vector<ScalarType, eigen> Vector<ScalarType, eigen>::sort() const
  {
      Vector<ScalarType, eigen> to_be_sorted = *this;
      std::sort( to_be_sorted.wrapped_.data(), to_be_sorted.wrapped_.data() + size() );
      return to_be_sorted;
  }

  //! Returns the indices of the first smallest coefficients
  template<class ScalarType>
  std::vector<size_t> Vector<ScalarType, eigen>::indices_of_smallest(size_t const number_of_indices)
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

  //! Returns the Euclidean norm
  template<class ScalarType> inline
  ScalarType Vector<ScalarType,eigen>::norm() const
  { return wrapped_.norm(); }

  template<class ScalarType> inline
  Vector<ScalarType,eigen> & Vector<ScalarType,eigen>::fill(ScalarType const& lambda)
  {
    wrapped_.fill(lambda);
    return *this;
  }

  template<class ScalarType> inline
  ScalarType Vector<ScalarType,eigen>::dot(Vector<ScalarType,eigen> const& v) const
  { return wrapped_.dot(v.wrapped_); }

  template<class ScalarType> inline
  ScalarType operator,(Vector<ScalarType,eigen> const& v, Vector<ScalarType, eigen> const & w)
  { return v.wrapped_.dot(w.wrapped_); }

  //======================
  // Operators
  //======================


  //! addition of another vector
  template<class ScalarType> inline
  Vector<ScalarType,eigen>& Vector<ScalarType,eigen>::operator+=(Vector<ScalarType,eigen> const& u)
  {
    wrapped_ += u.wrapped_;
    return *this;
  }

  //! substraction of another vector
  template<class ScalarType> inline
  Vector<ScalarType,eigen>& Vector<ScalarType,eigen>::operator-=(Vector<ScalarType,eigen> const& u)
  {
    wrapped_ -= u.wrapped_;
    return *this;
  }

  //! multiplication by a scalar
  template<class ScalarType> inline
  Vector<ScalarType,eigen>& Vector<ScalarType,eigen>::operator*=(ScalarType const& lambda)
  {
    wrapped_ *= lambda;
    return *this;
  }

  //! division by a scalar
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
  { return Vector<ScalarType,eigen>(*this) += u; }

  template<class ScalarType> inline
  Vector<ScalarType,eigen> Vector<ScalarType,eigen>::operator-(Vector<ScalarType,eigen> const& u) const
  { return Vector<ScalarType,eigen>(*this) -= u; }

  
  template<class ScalarType>
  Vector<double,eigen> piecewiseDivision(Vector<double,eigen> const& u, Vector<ScalarType,eigen> const& v)
  {
    if (u.size() != v.size())
      throw std::invalid_argument("Can only divide vectors of same size !");
    Vector<double,eigen> w(u.size());
    for (size_t i=0; i<u.size(); i++)
        w(i) = u(i) / v(i);
    return w;
  }

  //! Print a vector into a file
  template<class ScalarType>
  std::ostream & operator<<(std::ostream & fileToWrite, simol::Vector<ScalarType,eigen> const & vectorToRead)
  {
    for (size_t index = 0; index < vectorToRead.size(); ++index)
      fileToWrite << vectorToRead(index) << " ";
    return fileToWrite;
  }


}

#endif