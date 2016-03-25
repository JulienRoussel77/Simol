#ifndef SIMOL_SPARSEMATRIX_HPP
#define SIMOL_SPARSEMATRIX_HPP

#include "core/linalg/Vector.hpp"
#include "core/io/MatrixMarketFile.hpp"
//#include "quantchem/SparseTensor.hpp"
#include "eigen.hpp"

#include "dsaupd.hpp"
#include "dseupd.hpp"

#include <fstream>
#include <string>

namespace simol
{
    std::size_t
    getInd(std::size_t const M_disc,
           std::size_t const rowIndex,
           std::size_t const columnIndex);


  template<class ScalarType, template<class> class WrappedLibrary = eigen>
  class SparseMatrix;

  template<class ScalarType, template<class> class Library>
  std::ifstream & operator>>(std::ifstream & fileToRead,
                             SparseMatrix<ScalarType, Library> & matrixToWrite);

  template<class ScalarType>
  class SparseMatrix<ScalarType,eigen>
  {
    friend std::ifstream & operator>> <>(std::ifstream & fileToRead,
                                         SparseMatrix<ScalarType,eigen> & matrixToWrite);
    public:

      typedef typename eigen<ScalarType>::SparseMatrixType::InnerIterator iterator;

      operator typename eigen<ScalarType>::SparseMatrixType const () const
      { return wrapped_; }

      std::size_t numberOfRows() const;
      std::size_t numberOfColumns() const;
      ScalarType const operator()(std::size_t const rowIndex,
                                    std::size_t const columnIndex) const
     { return wrapped_.coeff(rowIndex, columnIndex); }

     ScalarType& operator()(std::size_t const rowIndex,
                                    std::size_t const columnIndex)
     { return wrapped_.coeffRef(rowIndex, columnIndex); }

      SparseMatrix<ScalarType,eigen>& operator+=(SparseMatrix<ScalarType,eigen> const& A);
      SparseMatrix<ScalarType,eigen>& operator-=(SparseMatrix<ScalarType,eigen> const& A);
      SparseMatrix<ScalarType,eigen>& operator*=(ScalarType const& scalar);
      SparseMatrix<ScalarType,eigen>& operator/=(ScalarType const& scalar);
      SparseMatrix<ScalarType,eigen> operator*(ScalarType const& scalar) const;
      SparseMatrix<ScalarType,eigen> operator/(ScalarType const& scalar) const;
      SparseMatrix<ScalarType,eigen> operator-() const;
      SparseMatrix<ScalarType,eigen> operator+(SparseMatrix<ScalarType,eigen> const& A) const;
      SparseMatrix<ScalarType,eigen> operator-(SparseMatrix<ScalarType,eigen> const& A) const;

      // TODO: write a non-pessimized version (return type may different in eigen)
      SparseMatrix adjoint() const
      { return SparseMatrix(wrapped_.adjoint()); }



    public:
      typename eigen<ScalarType>::SparseMatrixType const & wrapped() const
      { return wrapped_; }
    public:
      SparseMatrix(size_t const numberOfRows, size_t const numberOfColumns);
      SparseMatrix(MatrixMarketFile const & file);
      //SparseMatrix(const char * filename, std::size_t const size)
      //: SparseMatrix(filename, size) //wrapped_(size, size)
      //{
     	/*FILE* fichier = fopen(filename,"r" ); //ON ouvre le fichier en lecture seule
        typedef Eigen::Triplet<double,std::size_t> T;
	    std::vector<T> tripletList;
        std::ifstream in(filename); //Ouverture en mode lecture de "bdd.txt"
	    std::string ligne; //Création d'une chaine de caractere
	    int nbLignes = 0;

	    while(std::getline(in, ligne))
	    {
		    int i;
		    int j;
		    long double t;
		    int test = 0;
		    test = fscanf(fichier, "%d %d %Lf", &i , &j, &t);
		    assert(test>0);
		    double t0 =t;
		    tripletList.push_back(T(i,j,t0));

		    //On lit chaque ligne du fichier que l'on stoke dans "ligne"
	           nbLignes++;


	    }
	    in.close(); //On ferme le fichier
        wrapped_.setFromTriplets(tripletList.begin(), tripletList.end());*/
      //}
      SparseMatrix(std::string const & filename, std::size_t const size)
      : wrapped_(size, size)
      {
        FILE* fichier = fopen(filename.c_str(),"r" ); //ON ouvre le fichier en lecture seule
        typedef Eigen::Triplet<double,std::size_t> T;
        std::vector<T> tripletList;
          std::ifstream in(filename.c_str()); //Ouverture en mode lecture de "bdd.txt"
        std::string ligne; //Création d'une chaine de caractere
        int nbLignes = 0;

        while(std::getline(in, ligne))
        {
          int i;
          int j;
          long double t;
          int test = 0;
          test = fscanf(fichier, "%d %d %Lf", &i , &j, &t);
          assert(test>0);
          double t0 =t;
          tripletList.push_back(T(i,j,t0));

          //On lit chaque ligne du fichier que l'on stoke dans "ligne"
              nbLignes++;
        }
        in.close(); //On ferme le fichier
          wrapped_.setFromTriplets(tripletList.begin(), tripletList.end());
      }

      SparseMatrix(Vector<double, eigen> u, std::size_t numberOfRows, std::size_t numberOfColumns);
      SparseMatrix(const typename eigen<ScalarType>::SparseMatrixType::AdjointReturnType& A);
    public:
      typename eigen<ScalarType>::SparseMatrixType wrapped_;
  };




  template<class ScalarType> inline
  SparseMatrix<ScalarType,eigen>::SparseMatrix(size_t const numberOfRows,
                                               size_t const numberOfColumns)
  : wrapped_(numberOfRows,numberOfColumns)
  {}

  template<class ScalarType> inline
  SparseMatrix<ScalarType,eigen>::SparseMatrix(MatrixMarketFile const & file)
  : wrapped_(file.numberOfRows(), file.numberOfColumns())
  {
    std::vector< Eigen::Triplet<ScalarType, std::size_t> > nonzeros(file.numberOfNonzeros());
    for(size_t nonzeroIndex = 0; nonzeroIndex < file.numberOfNonzeros(); ++nonzeroIndex)
    {
        int rowIndex;
        int columnIndex;
        ScalarType nonzero;
        fscanf(file.content(), "%d %d %lg\n", &rowIndex, &columnIndex, &nonzero);
        nonzeros.push_back(Eigen::Triplet<ScalarType,size_t>(rowIndex, columnIndex, nonzero));
    }
    wrapped_.setFromTriplets(nonzeros.begin(),nonzeros.end());
  }

  template<class ScalarType>
  SparseMatrix<ScalarType, eigen>::SparseMatrix(Vector<double, eigen> u, std::size_t numberOfRows, std::size_t numberOfColumns):
    SparseMatrix(numberOfRows, numberOfColumns)
  {

    //wrapped_.reshape(u, numberOfRows, numberOfColumns);

    if (u.size() != numberOfColumns * numberOfRows)
      throw std::invalid_argument("Matrix reshape : sizes incoherent !");
    for (int j=0; j<numberOfColumns; j++)
      for (SparseMatrix<double>::iterator it(*this,j); it; ++it)
        wrapped_.coeffRef(it.row(), j) = u(it.row()*u.size() + j);
  }

  template<class ScalarType>
  SparseMatrix<ScalarType, eigen>::SparseMatrix(const typename eigen<ScalarType>::SparseMatrixType::AdjointReturnType& A)
  :wrapped_(A)
  {}

  template<class ScalarType> inline std::size_t
  SparseMatrix<ScalarType,eigen>::numberOfRows() const
  { return wrapped_.rows(); }

  template<class ScalarType> inline std::size_t
  SparseMatrix<ScalarType,eigen>::numberOfColumns() const
  { return wrapped_.cols(); }


  /*template<class ScalarType>
  SparseMatrix<ScalarType,eigen> SparseMatrix<ScalarType,eigen>::operator*=(ScalarType const& scalar)
  {
    return SparseMatrix (wrapped_ * scalar);
  }

  template<class ScalarType>
  SparseMatrix<ScalarType,eigen> SparseMatrix<ScalarType,eigen>::operator*=(ScalarType const& scalar)
  {
    return SparseMatrix (wrapped_ * scalar);
  }*/





  template<class ScalarType, template<class> class Library>
  std::ifstream & operator>>(std::ifstream & fileToRead, SparseMatrix<ScalarType,Library> & matrixToWrite)
  {
    std::vector< Eigen::Triplet<ScalarType,size_t> > nonzeros;

    size_t rowIndex, columnIndex;
    ScalarType value;

    while(fileToRead >> rowIndex >> columnIndex >> value)
      nonzeros.push_back(Eigen::Triplet<ScalarType,size_t>(rowIndex,columnIndex,value));

    matrixToWrite.wrapped_.setFromTriplets(nonzeros.begin(),nonzeros.end());

    return fileToRead;
  }

  template<class ScalarType>
  Vector<ScalarType, eigen> operator*(SparseMatrix<ScalarType, eigen> const & matrix,
                                      Vector<ScalarType, eigen> const & vector)
  {
      Vector<ScalarType, eigen> prod(matrix.numberOfRows());
      prod.wrapped_ = matrix.wrapped_.template selfadjointView<Eigen::Upper>() * vector.wrapped_;
      return prod;
  }

  /*template<class ScalarType>
  SparseMatrix<ScalarType, eigen> operator+(SparseMatrix<ScalarType, eigen> const & A,
                                      SparseMatrix<ScalarType, eigen> const & B)
  {
      SparseMatrix<ScalarType, eigen> C(A.numberOfRows(), A.numberOfColumns());
      C.wrapped_ = A.wrapped_ - B.wrapped_;
      return C;
  }*/

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
  SparseMatrix<ScalarType, eigen> speye(size_t const nbOfRows, size_t const nbOfColumns)
  {
    SparseMatrix<ScalarType, eigen> A(nbOfRows, nbOfColumns);
    for (int i=0; i < std::min(nbOfRows, nbOfColumns); i++)
      A(i,i) = 1;
    //A.wrapped_.setIdentity();
    return A;
  }

  template<class ScalarType>
  SparseMatrix<ScalarType, eigen> spzero(size_t const nbOfRows, size_t const nbOfColumns)
  {
    SparseMatrix<ScalarType, eigen> A(nbOfRows, nbOfColumns);
    return A;
  }

  template<class ScalarType>
  SparseMatrix<ScalarType,eigen> operator*(ScalarType const& scalar, SparseMatrix<ScalarType,eigen> const& A)
  { return SparseMatrix<ScalarType,eigen>(typename eigen<ScalarType>::SparseMatrixType(scalar*A.wrapped_)); }

}


#endif
