#ifndef SIMOL_SPARSEMATRIX_EIGEN_HPP
#define SIMOL_SPARSEMATRIX_EIGEN_HPP

#include "DenseMatrix.hpp"

namespace simol
{
  template<class ScalarType>
  class SparseMatrix<ScalarType, eigen>
  {
      friend std::ifstream & operator>> <>(std::ifstream & fileToRead, SparseMatrix<ScalarType, eigen> & matrixToWrite);

    public:
      typedef typename eigen<ScalarType>::SparseMatrixType::InnerIterator iterator;

      explicit SparseMatrix(size_t const numberOfRows, size_t const numberOfColumns);
      explicit SparseMatrix(MatrixMarketFile const & file);
      explicit SparseMatrix(std::string const & filename, std::size_t const size);
      explicit SparseMatrix(Vector<double, eigen> u, std::size_t numberOfRows, std::size_t numberOfColumns);
      SparseMatrix(const typename eigen<ScalarType>::SparseMatrixType::AdjointReturnType& A);

      operator typename eigen<ScalarType>::SparseMatrixType const & () const;

      std::size_t numberOfRows() const;
      std::size_t numberOfColumns() const;
      std::size_t nonZeros() const;

      ScalarType const operator()(std::size_t const rowIndex, std::size_t const columnIndex) const;
      ScalarType & operator()(std::size_t const rowIndex, std::size_t const columnIndex);
      ScalarType& insert(std::size_t const rowIndex, std::size_t const columnIndex);

      SparseMatrix<ScalarType, eigen>& operator+=(SparseMatrix<ScalarType, eigen> const& A);
      SparseMatrix<ScalarType, eigen>& operator-=(SparseMatrix<ScalarType, eigen> const& A);
      SparseMatrix<ScalarType, eigen>& operator*=(ScalarType const& scalar);
      SparseMatrix<ScalarType, eigen>& operator/=(ScalarType const& scalar);
      SparseMatrix<ScalarType, eigen> operator*(ScalarType const& scalar) const;
      SparseMatrix<ScalarType, eigen> operator/(ScalarType const& scalar) const;
      SparseMatrix<ScalarType, eigen> operator-() const;
      SparseMatrix<ScalarType, eigen> operator+(SparseMatrix<ScalarType, eigen> const& A) const;
      SparseMatrix<ScalarType, eigen> operator-(SparseMatrix<ScalarType, eigen> const& A) const;

      SparseMatrix adjoint() const;

      typename eigen<ScalarType>::SparseMatrixType const & wrapped() const;

      virtual DenseMatrix<ScalarType, eigen> dense() const; 

    public:
      typename eigen<ScalarType>::SparseMatrixType wrapped_;
  };

  //! Conversion to wrapped matrix
  template<class ScalarType> inline
  SparseMatrix<ScalarType, eigen>::operator typename eigen<ScalarType>::SparseMatrixType const & () const
  { return wrapped_; }

  //! Returns coefficient
  template<class ScalarType> inline
  ScalarType const SparseMatrix<ScalarType, eigen>::operator()(std::size_t const rowIndex, std::size_t const columnIndex) const
  { return wrapped_.coeff(rowIndex, columnIndex); }

  //! Modify coefficient
  template<class ScalarType> inline
  ScalarType& SparseMatrix<ScalarType, eigen>::operator()(std::size_t const rowIndex, std::size_t const columnIndex)
  { return wrapped_.coeffRef(rowIndex, columnIndex); }

  //! Create coefficient
  template<class ScalarType> inline
  ScalarType& SparseMatrix<ScalarType, eigen>::insert(std::size_t const rowIndex, std::size_t const columnIndex)
  { return wrapped_.insert(rowIndex, columnIndex); }

  //! Returns adjoint matrix
  template<class ScalarType> inline
  typename eigen<ScalarType>::SparseMatrixType const & SparseMatrix<ScalarType, eigen>::wrapped() const
  { return wrapped_; }

  //! Returns adjoint matrix
  // TODO: write a non-pessimized version (return type may different in eigen)
  template<class ScalarType> inline
  SparseMatrix<ScalarType, eigen> SparseMatrix<ScalarType, eigen>::adjoint() const
  { return SparseMatrix(wrapped_.adjoint()); }


    //! Returns a dense matrix equal to the sparse matrix
    template<typename ScalarType> inline
    DenseMatrix<ScalarType, eigen> SparseMatrix<ScalarType, eigen>::dense() const
    {
       DenseMatrix<ScalarType, eigen> M(wrapped_.rows(), wrapped_.cols());
       M.wrapped_ = wrapped_.template triangularView<Eigen::Upper>();
       M.wrapped_ += wrapped_.template triangularView<Eigen::StrictlyLower>();
        return M;
    }





  //! Construction from a file
  template<class ScalarType>
  SparseMatrix<ScalarType, eigen>::SparseMatrix(std::string const & filename, std::size_t const size)
    : wrapped_(size, size)
  {
    FILE* fichier = fopen(filename.c_str(), "r" ); //ON ouvre le fichier en lecture seule
    typedef Eigen::Triplet<double, std::size_t> T;
    std::vector<T> tripletList;
    std::ifstream in(filename.c_str()); //Ouverture en mode lecture de "bdd.txt"
    std::string ligne; //Création d'une chaine de caractere
    int nbLignes = 0;

    while(std::getline(in, ligne))
    {
      int i;
      int j;
      long double t;
      fscanf(fichier, "%d %d %Lf", &i , &j, &t);
      double t0 = t;
      tripletList.push_back(T(i, j, t0));

      //On lit chaque ligne du fichier que l'on stoke dans "ligne"
      nbLignes++;
    }
    in.close(); //On ferme le fichier
    wrapped_.setFromTriplets(tripletList.begin(), tripletList.end());
  }


  //! Construction from a given size
  template<class ScalarType> inline
  SparseMatrix<ScalarType, eigen>::SparseMatrix(size_t const numberOfRows, size_t const numberOfColumns)
    : wrapped_(numberOfRows, numberOfColumns)
  {}

  //! Construction from a Matrix Market file
  template<class ScalarType> inline
  SparseMatrix<ScalarType, eigen>::SparseMatrix(MatrixMarketFile const & file)
    : wrapped_(file.numberOfRows(), file.numberOfColumns())
  {
    std::vector< Eigen::Triplet<ScalarType, std::size_t> > nonzeros(file.numberOfNonzeros());
    for(size_t nonzeroIndex = 0; nonzeroIndex < file.numberOfNonzeros(); ++nonzeroIndex)
    {
      int rowIndex;
      int columnIndex;
      ScalarType nonzero;
      fscanf(file.content(), "%d %d %lg\n", &rowIndex, &columnIndex, &nonzero);
      nonzeros.push_back(Eigen::Triplet<ScalarType, size_t>(rowIndex, columnIndex, nonzero));
    }
    wrapped_.setFromTriplets(nonzeros.begin(), nonzeros.end());
  }

  //! Construction from a vector
  template<class ScalarType>
  SparseMatrix<ScalarType, eigen>::SparseMatrix(Vector<double, eigen> u, std::size_t numberOfRows, std::size_t numberOfColumns):
    SparseMatrix(numberOfRows, numberOfColumns)
  {
    if (u.size() != numberOfColumns * numberOfRows)
      throw std::invalid_argument("Matrix reshape : sizes incoherent !");
    for (std::size_t j = 0; j < numberOfColumns; j++)
      for (SparseMatrix<double>::iterator it(*this, j); it; ++it)
        wrapped_.coeffRef(it.row(), j) = u(it.row() * u.size() + j);
  }

  template<class ScalarType>
  SparseMatrix<ScalarType, eigen>::SparseMatrix(const typename eigen<ScalarType>::SparseMatrixType::AdjointReturnType& A)
    : wrapped_(A)
  {}

  template<class ScalarType> inline std::size_t
  SparseMatrix<ScalarType, eigen>::numberOfRows() const
  { return wrapped_.rows(); }

  template<class ScalarType> inline std::size_t
  SparseMatrix<ScalarType, eigen>::numberOfColumns() const
  { return wrapped_.cols(); }

  template<class ScalarType> inline std::size_t
  SparseMatrix<ScalarType, eigen>::nonZeros() const
  { return wrapped_.nonZeros(); }


  template<class ScalarType, template<class> class Library>
  std::ifstream & operator>>(std::ifstream & fileToRead, SparseMatrix<ScalarType, Library> & matrixToWrite)
  {
    std::vector< Eigen::Triplet<ScalarType, size_t> > nonzeros;

    size_t rowIndex, columnIndex;
    ScalarType value;

    while(fileToRead >> rowIndex >> columnIndex >> value)
      nonzeros.push_back(Eigen::Triplet<ScalarType, size_t>(rowIndex, columnIndex, value));

    matrixToWrite.wrapped_.setFromTriplets(nonzeros.begin(), nonzeros.end());

    return fileToRead;
  }

  // Does a matrix by vector product
  template<class ScalarType>
  Vector<ScalarType, eigen> operator*(SparseMatrix<ScalarType, eigen> matrix, Vector<ScalarType, eigen> const & vector)
  {
    Vector<ScalarType, eigen> prod(matrix.wrapped_.rows());
    prod.wrapped_ = matrix.wrapped_ * vector.wrapped_;
    return prod;
  }


  //======================
  // Operators
  //======================


  template<class ScalarType> inline
  SparseMatrix<ScalarType>& SparseMatrix<ScalarType>::operator+=(SparseMatrix<ScalarType> const& A)
  {
    wrapped_ += A.wrapped_;
    return *this;
  }

  template<class ScalarType> inline
  SparseMatrix<ScalarType>& SparseMatrix<ScalarType>::operator-=(SparseMatrix<ScalarType> const& A)
  {
    wrapped_ -= A.wrapped_;
    return *this;
  }

  template<class ScalarType> inline
  SparseMatrix<ScalarType>& SparseMatrix<ScalarType>::operator*=(ScalarType const& scalar)
  {
    wrapped_ *= scalar;
    return *this;
  }

  template<class ScalarType> inline
  SparseMatrix<ScalarType>& SparseMatrix<ScalarType>::operator/=(ScalarType const& scalar)
  {
    wrapped_ /= scalar;
    return *this;
  }

  template<class ScalarType> inline
  SparseMatrix<ScalarType> SparseMatrix<ScalarType>::operator*(ScalarType const& scalar) const
  {
    Vector<ScalarType> A(*this);
    A.wrapped_ *= scalar;
    return A;
  }

  template<class ScalarType> inline
  SparseMatrix<ScalarType> SparseMatrix<ScalarType>::operator/(ScalarType const& scalar) const
  {
    Vector<ScalarType> A(*this);
    A.wrapped_ /= scalar;
    return A;
  }

  template<class ScalarType> inline
  SparseMatrix<ScalarType> SparseMatrix<ScalarType>::operator-() const
  {
    Vector<ScalarType> A(*this);
    A.wrapped_ *= -1;
    return A;
  }

  template<class ScalarType> inline
  SparseMatrix<ScalarType> SparseMatrix<ScalarType>::operator+(SparseMatrix<ScalarType> const& A) const
  { return SparseMatrix<ScalarType>(*this) += A; }

  template<class ScalarType> inline
  SparseMatrix<ScalarType> SparseMatrix<ScalarType>::operator-(SparseMatrix<ScalarType> const& A) const
  { return SparseMatrix<ScalarType>(*this) -= A; }


  template<class ScalarType>
  SparseMatrix<ScalarType> speye(size_t const nbOfRows, size_t const nbOfColumns)
  {
    SparseMatrix<ScalarType> A(nbOfRows, nbOfColumns);
    for (std::size_t i = 0; i < std::min(nbOfRows, nbOfColumns); i++)
      A(i, i) = 1;
    //A.insert(i,i) = 1;
    return A;
  }

  template<class ScalarType>
  SparseMatrix<ScalarType, eigen> spzero(size_t const nbOfRows, size_t const nbOfColumns)
  {
    SparseMatrix<ScalarType, eigen> A(nbOfRows, nbOfColumns);
    return A;
  }

  template<class ScalarType>
  SparseMatrix<ScalarType, eigen> operator*(ScalarType const& scalar, SparseMatrix<ScalarType, eigen> const& A)
  { return SparseMatrix<ScalarType, eigen>(typename eigen<ScalarType>::SparseMatrixType(scalar * A.wrapped_)); }

}

#endif
