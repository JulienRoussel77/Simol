#ifndef SIMOL_SPARSEMATRIX_EIGEN_HPP
#define SIMOL_SPARSEMATRIX_EIGEN_HPP

namespace simol
{
  template<class Scalar>
  class SparseMatrix<Scalar, eigen>
  {
      friend std::ifstream & operator>> <>(std::ifstream & fileToRead, SparseMatrix<Scalar, eigen> & matrixToWrite);

    public:


      operator typename eigen::SparseMatrix<Scalar> const & () const;

      explicit SparseMatrix(size_t const numberOfRows, size_t const numberOfColumns);
      explicit SparseMatrix(MatrixMarketFile const & file);
      explicit SparseMatrix(std::string const & filename, std::size_t const size);
      explicit SparseMatrix(Vector<double, eigen> u, std::size_t numberOfRows, std::size_t numberOfColumns);
      SparseMatrix(const typename eigen::SparseMatrix<Scalar>::AdjointReturnType& A);

      std::size_t number_of_rows() const;
      std::size_t number_of_columns() const;

      Scalar const operator()(std::size_t const rowIndex, std::size_t const columnIndex) const;
      Scalar & operator()(std::size_t const rowIndex, std::size_t const columnIndex);
      Scalar & insert(std::size_t const rowIndex, std::size_t const columnIndex);

      SparseMatrix adjoint() const;

      SparseMatrix<Scalar, eigen>& operator+=(SparseMatrix<Scalar, eigen> const& A);
      SparseMatrix<Scalar, eigen>& operator-=(SparseMatrix<Scalar, eigen> const& A);
      SparseMatrix<Scalar, eigen>& operator*=(Scalar const& scalar);
      SparseMatrix<Scalar, eigen>& operator/=(Scalar const& scalar);
      SparseMatrix<Scalar, eigen> operator*(Scalar const& scalar) const;
      SparseMatrix<Scalar, eigen> operator/(Scalar const& scalar) const;
      SparseMatrix<Scalar, eigen> operator-() const;
      SparseMatrix<Scalar, eigen> operator+(SparseMatrix<Scalar, eigen> const& A) const;
      SparseMatrix<Scalar, eigen> operator-(SparseMatrix<Scalar, eigen> const& A) const;


      typename eigen::SparseMatrix<Scalar> const & wrapped() const;

    public:
      typedef typename eigen::SparseMatrix<Scalar>::InnerIterator iterator;
      typename eigen::SparseMatrix<Scalar> wrapped_;
  };

  //! Conversion to wrapped matrix
  template<class Scalar> inline
  SparseMatrix<Scalar, eigen>::operator typename eigen::SparseMatrix<Scalar> const & () const
  { return wrapped_; }

  //! Returns coefficient
  template<class Scalar> inline
  Scalar const SparseMatrix<Scalar, eigen>::operator()(std::size_t const rowIndex, std::size_t const columnIndex) const
  { return wrapped_.coeff(rowIndex, columnIndex); }

  //! Modify coefficient
  template<class Scalar> inline
  Scalar& SparseMatrix<Scalar, eigen>::operator()(std::size_t const rowIndex, std::size_t const columnIndex)
  { return wrapped_.coeffRef(rowIndex, columnIndex); }

  //! Create coefficient
  template<class Scalar> inline
  Scalar& SparseMatrix<Scalar, eigen>::insert(std::size_t const rowIndex, std::size_t const columnIndex)
  { return wrapped_.insert(rowIndex, columnIndex); }

  //! Convert into wrapped matrix type
  template<class Scalar> inline
  typename eigen::SparseMatrix<Scalar> const & SparseMatrix<Scalar, eigen>::wrapped() const
  { return wrapped_; }

  //! Returns adjoint matrix
  // TODO: write a non-pessimized version (return type may different in eigen)
  template<typename Scalar> inline
  SparseMatrix<Scalar, eigen> SparseMatrix<Scalar, eigen>::adjoint() const
  { return SparseMatrix(eigen::adjoint(wrapped_)); }


  //! Construction from a file
  template<class Scalar>
  SparseMatrix<Scalar, eigen>::SparseMatrix(std::string const & filename, std::size_t const size)
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
  template<class Scalar> inline
  SparseMatrix<Scalar, eigen>::SparseMatrix(size_t const numberOfRows, size_t const numberOfColumns)
    : wrapped_(numberOfRows, numberOfColumns)
  {}

  //! Construction from a Matrix Market file
  template<class Scalar> inline
  SparseMatrix<Scalar, eigen>::SparseMatrix(MatrixMarketFile const & file)
    : wrapped_(file.numberOfRows(), file.numberOfColumns())
  {
    std::vector< Eigen::Triplet<Scalar, std::size_t> > nonzeros(file.numberOfNonzeros());
    for(size_t nonzeroIndex = 0; nonzeroIndex < file.numberOfNonzeros(); ++nonzeroIndex)
    {
      int rowIndex;
      int columnIndex;
      Scalar nonzero;
      fscanf(file.content(), "%d %d %lg\n", &rowIndex, &columnIndex, &nonzero);
      nonzeros.push_back(Eigen::Triplet<Scalar, size_t>(rowIndex, columnIndex, nonzero));
    }
    wrapped_.setFromTriplets(nonzeros.begin(), nonzeros.end());
  }

  //! Construction from a vector
  template<class Scalar>
  SparseMatrix<Scalar, eigen>::SparseMatrix(Vector<double, eigen> u, std::size_t numberOfRows, std::size_t numberOfColumns):
    SparseMatrix(numberOfRows, numberOfColumns)
  {
    if (u.size() != numberOfColumns * numberOfRows)
      throw std::invalid_argument("Matrix reshape : sizes incoherent !");
    for (std::size_t j = 0; j < numberOfColumns; j++)
      for (SparseMatrix<double>::iterator it(*this, j); it; ++it)
        wrapped_.coeffRef(it.row(), j) = u(it.row() * u.size() + j);
  }

  template<class Scalar>
  SparseMatrix<Scalar, eigen>::SparseMatrix(const typename eigen::SparseMatrix<Scalar>::AdjointReturnType& A)
    : wrapped_(A)
  {}

  //! Returns the number of rows
  template<class Scalar> inline
  std::size_t SparseMatrix<Scalar, eigen>::number_of_rows() const
  { return eigen::number_of_rows(wrapped_); }

  //! Returns the number of columns
  template<class Scalar> inline
  std::size_t SparseMatrix<Scalar, eigen>::number_of_columns() const
  { return eigen::number_of_columns(wrapped_); }



  template<typename Scalar>
  std::ifstream & operator>>(std::ifstream & fileToRead, SparseMatrix<Scalar, eigen> & matrixToWrite)
  {
    std::vector< Eigen::Triplet<Scalar, size_t> > nonzeros;

    size_t rowIndex, columnIndex;
    Scalar value;

    while(fileToRead >> rowIndex >> columnIndex >> value)
      nonzeros.push_back(Eigen::Triplet<Scalar, size_t>(rowIndex, columnIndex, value));

    matrixToWrite.wrapped_.setFromTriplets(nonzeros.begin(), nonzeros.end());

    return fileToRead;
  }

  // TODO: attention version symetrique
  template<class Scalar>
  Vector<Scalar, eigen> operator*(SparseMatrix<Scalar, eigen> const & matrix,
                                  Vector<Scalar, eigen> const & vector)
  {
    Vector<Scalar, eigen> prod(matrix.number_of_rows());
    prod.wrapped_ = matrix.wrapped_.template selfadjointView<Eigen::Upper>() * vector.wrapped_;
    return prod;
  }


  //======================
  // Operators
  //======================


  template<class Scalar> inline
  SparseMatrix<Scalar>& SparseMatrix<Scalar>::operator+=(SparseMatrix<Scalar> const& A)
  {
    wrapped_ += A.wrapped_;
    return *this;
  }

  template<class Scalar> inline
  SparseMatrix<Scalar>& SparseMatrix<Scalar>::operator-=(SparseMatrix<Scalar> const& A)
  {
    wrapped_ -= A.wrapped_;
    return *this;
  }

  template<class Scalar> inline
  SparseMatrix<Scalar>& SparseMatrix<Scalar>::operator*=(Scalar const& scalar)
  {
    wrapped_ *= scalar;
    return *this;
  }

  template<class Scalar> inline
  SparseMatrix<Scalar>& SparseMatrix<Scalar>::operator/=(Scalar const& scalar)
  {
    wrapped_ /= scalar;
    return *this;
  }

  template<class Scalar> inline
  SparseMatrix<Scalar> SparseMatrix<Scalar>::operator*(Scalar const& scalar) const
  {
    Vector<Scalar> A(*this);
    A.wrapped_ *= scalar;
    return A;
  }

  template<class Scalar> inline
  SparseMatrix<Scalar> SparseMatrix<Scalar>::operator/(Scalar const& scalar) const
  {
    Vector<Scalar> A(*this);
    A.wrapped_ /= scalar;
    return A;
  }

  template<class Scalar> inline
  SparseMatrix<Scalar> SparseMatrix<Scalar>::operator-() const
  {
    Vector<Scalar> A(*this);
    A.wrapped_ *= -1;
    return A;
  }

  template<class Scalar> inline
  SparseMatrix<Scalar> SparseMatrix<Scalar>::operator+(SparseMatrix<Scalar> const& A) const
  { return SparseMatrix<Scalar>(*this) += A; }

  template<class Scalar> inline
  SparseMatrix<Scalar> SparseMatrix<Scalar>::operator-(SparseMatrix<Scalar> const& A) const
  { return SparseMatrix<Scalar>(*this) -= A; }


  template<class Scalar>
  SparseMatrix<Scalar> speye(size_t const nbOfRows, size_t const nbOfColumns)
  {
    SparseMatrix<Scalar> A(nbOfRows, nbOfColumns);
    for (std::size_t i = 0; i < std::min(nbOfRows, nbOfColumns); i++)
      A(i, i) = 1;
    //A.insert(i,i) = 1;
    return A;
  }

  template<class Scalar>
  SparseMatrix<Scalar, eigen> spzero(size_t const nbOfRows, size_t const nbOfColumns)
  {
    SparseMatrix<Scalar, eigen> A(nbOfRows, nbOfColumns);
    return A;
  }

  template<class Scalar>
  SparseMatrix<Scalar, eigen> operator*(Scalar const& scalar, SparseMatrix<Scalar, eigen> const& A)
  { return SparseMatrix<Scalar, eigen>(typename eigen::SparseMatrix<Scalar>(scalar * A.wrapped_)); }

}

#endif
