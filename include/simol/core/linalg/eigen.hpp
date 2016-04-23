#ifndef EIGEN_HPP
#define EIGEN_HPP

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wsign-compare"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCore>
#include <Eigen/SVD>
#pragma GCC diagnostic pop


namespace simol
{

  struct eigen
  {

    //--- aliases ---

    template<typename Scalar> using DenseMatrix      = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    template<typename Scalar> using DenseAdjointd    = typename DenseMatrix<Scalar>::AdjointReturnType;
    template<typename Scalar> using DenseColumn      = typename DenseMatrix<Scalar>::ColXpr;
    template<typename Scalar> using DenseBlock       = Eigen::Block<DenseMatrix<Scalar>>;
    template<typename Scalar> using DenseBlock_const = Eigen::Block<DenseMatrix<Scalar> const>;
    template<typename Scalar> using Vector           = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    template<typename Scalar> using VectorMap        = Eigen::Map<Vector<Scalar>>;
    template<typename Scalar> using SparseMatrix     = Eigen::SparseMatrix<Scalar>;
    template<typename Scalar> using SVD              = Eigen::JacobiSVD<DenseMatrix<Scalar>>;

    //--- static functions ---

    template<typename Scalar> static void multiply(Vector<Scalar> & vector, Scalar const scalar);
    template<typename Scalar> static void add(Vector<Scalar> & lhs, Vector<Scalar> const & rhs);
    template<typename Scalar> static void fill(Vector<Scalar> & vector, Scalar const scalar);
    
    template<typename Scalar> static std::size_t size(Vector<Scalar> const & vector);
    template<typename Scalar> static std::size_t min_index(Vector<Scalar> const & vector);
    template<typename Scalar> static std::size_t number_of_rows(DenseMatrix<Scalar> const & matrix);
    template<typename Scalar> static std::size_t number_of_rows(SparseMatrix<Scalar> const & matrix);
    template<typename Scalar> static std::size_t number_of_columns(DenseMatrix<Scalar> const & matrix);
    template<typename Scalar> static std::size_t number_of_columns(SparseMatrix<Scalar> const & matrix);

    template<typename Scalar> static Scalar min(Vector<Scalar> const & vector);
    template<typename Scalar> static Scalar max(Vector<Scalar> const & vector);
    template<typename Scalar> static Scalar norm(Vector<Scalar> const & vector);
    template<typename Scalar> static Scalar inner_product(Vector<Scalar> const & lhs, Vector<Scalar> const & rhs);
    template<typename Scalar> static Scalar determinant(DenseMatrix<Scalar> const & matrix);
    template<typename Scalar> static Scalar trace(DenseMatrix<Scalar> const & matrix);

    template<typename Scalar> static Vector<Scalar> sort(Vector<Scalar> const & vector);
    template<typename Scalar> static Vector<Scalar> subvector(Vector<Scalar> const & vector, std::size_t start, std::size_t length);
    template<typename Scalar> static Vector<Scalar> Zero(std::size_t const length);

    template<typename Scalar> static DenseMatrix<Scalar> Zero(std::size_t const numberOfRows, std::size_t const numberOfColumns);
    template<typename Scalar> static DenseMatrix<Scalar> Identity(std::size_t const dimension);
    template<typename Scalar> static DenseMatrix<Scalar> adjoint(DenseMatrix<Scalar> const & matrix);
    template<typename Scalar> static DenseMatrix<Scalar> inverse(DenseMatrix<Scalar> const & matrix);

    template<typename Scalar> static DenseBlock<Scalar> block(DenseMatrix<Scalar> & matrix, std::size_t const startRow, std::size_t const startCol, std::size_t const blockRows, std::size_t const blockCols);
    template<typename Scalar> static DenseBlock_const<Scalar> block(std::size_t const startRow, std::size_t const startCol, std::size_t const blockRows, std::size_t const blockCols);
   
   template<typename Scalar> static SparseMatrix<Scalar> adjoint(SparseMatrix<Scalar> const & matrix);

 };

  //! Multiply a vector by a scalar
  template<typename Scalar> inline
  void eigen::multiply(eigen::Vector<Scalar> & vector, 
                      Scalar const scalar)
  { vector *= scalar; }
  
  //! Add a vector to another one
  template<typename Scalar> inline
  void eigen::add(eigen::Vector<Scalar> & lhs, 
                  eigen::Vector<Scalar> const & rhs)
  { lhs += rhs; }
  
  //! Fill a vector with a scalar
  template<typename Scalar> inline
  void eigen::fill(eigen::Vector<Scalar> & vector, 
                   Scalar const scalar)
  { vector.fill(scalar); }
  
  //! Returns the size of a vector
  template<typename Scalar> inline
  std::size_t eigen::size(eigen::Vector<Scalar> const & vector)
  { return vector.size(); }

  //! Returns the index of the minimum coefficient of a vector
  template<typename Scalar> inline
  std::size_t eigen::min_index(eigen::Vector<Scalar> const & vector)
  {
    std::size_t index;
    vector.minCoeff(&index);
    return index;
  }

  //! Returns the number of rows of a dense matrix
  template<typename Scalar> inline
  std::size_t eigen::number_of_rows(eigen::DenseMatrix<Scalar> const & matrix)
  { return matrix.rows(); }

  //! Returns the number of rows of a sparse matrix
  template<typename Scalar> inline
  std::size_t eigen::number_of_rows(eigen::SparseMatrix<Scalar> const & matrix)
  { return matrix.rows(); }

  //! Returns the number of columns of a dense matrix
  template<typename Scalar> inline
  std::size_t eigen::number_of_columns(eigen::DenseMatrix<Scalar> const & matrix)
  { return matrix.cols(); }

  //! Returns the number of columns of a sparse matrix
  template<typename Scalar> inline
  std::size_t eigen::number_of_columns(eigen::SparseMatrix<Scalar> const & matrix)
  { return matrix.cols(); }

  //! Returns the minimum coefficient of a vector
  template<typename Scalar> inline
  Scalar eigen::min(eigen::Vector<Scalar> const & vector)
  { return vector.minCoeff(); }

  //! Returns the maximum coefficient of a vector
  template<typename Scalar> inline
  Scalar eigen::max(eigen::Vector<Scalar> const & vector)
  { return vector.maxCoeff(); }

  //! Returns the norm of a vector
  template<typename Scalar> inline
  Scalar eigen::norm(eigen::Vector<Scalar> const & vector)
  { return vector.norm(); }

  //! Return the trace of a dense matrix
  template<typename Scalar> inline
  Scalar eigen::trace(eigen::DenseMatrix<Scalar> const & matrix)
  { return matrix.trace(); }

  //! Returns the determinant of a dense matrix
  template<typename Scalar> inline
  Scalar eigen::determinant(eigen::DenseMatrix<Scalar> const & matrix)
  { return matrix.determinant(); }


  //! Returns a sorted copy of a vector
  template<typename Scalar> inline
  typename eigen::Vector<Scalar> eigen::sort(eigen::Vector<Scalar> const & vector)
  {
    typename eigen::Vector<Scalar> to_be_sorted = vector;
    std::sort( to_be_sorted.data(), to_be_sorted.data() + to_be_sorted.size() );
    return to_be_sorted;
  }

  //! Returns the null vector
  template<typename Scalar> inline
  typename eigen::Vector<Scalar> eigen::Zero(std::size_t const length)
  { return eigen::Vector<Scalar>(eigen::Vector<Scalar>::Zero(length)); }

  // TODO: could be optimized by calling a suitable function of Eigen
  //! Returns a subvector of a vector
  template<typename Scalar>
  typename eigen::Vector<Scalar> eigen::subvector(eigen::Vector<Scalar> const & vector,
      std::size_t start,
      std::size_t length)
  {
    eigen::Vector<Scalar> subvec(length);
    for (std::size_t index = 0; index < length; ++index)
      subvec(index) = vector(index + start);
    return subvec;
  }

  //! Returns the inner product between two vectors
  template<typename Scalar> inline
  Scalar eigen::inner_product(eigen::Vector<Scalar> const & lhs,
                              eigen::Vector<Scalar> const & rhs)
  { return lhs.dot(rhs); }

  //! Returns a null dense matrix
  template<typename Scalar> inline
  typename eigen::DenseMatrix<Scalar> eigen::Zero(std::size_t const numberOfRows, size_t const numberOfColumns)
  { return eigen::DenseMatrix<Scalar>(eigen::DenseMatrix<Scalar>::Zero(numberOfRows, numberOfColumns)); }

  //! Returns an identity dense matrix
  template<typename Scalar> inline
  typename eigen::DenseMatrix<Scalar> eigen::Identity(std::size_t const dimension)
  { return eigen::DenseMatrix<Scalar>(eigen::DenseMatrix<Scalar>::Identity(dimension, dimension)); }

  // TODO: write a non-pessimized version (return type may different in eigen)
  //! Adjoint of a dense matrix
  template<typename Scalar> inline
  typename eigen::DenseMatrix<Scalar> eigen::adjoint(eigen::DenseMatrix<Scalar> const & matrix)
  { return eigen::DenseMatrix<Scalar>(matrix.adjoint()); }
  
  //! Returns the inverse of a dense matrix
  template<typename Scalar> inline
  typename eigen::DenseMatrix<Scalar> eigen::inverse(eigen::DenseMatrix<Scalar> const & matrix)
  { return eigen::DenseMatrix<Scalar>(matrix.inverse()); }

  //! Returns a mutable submatrix
  template<typename Scalar> inline
  typename eigen::DenseBlock<Scalar> block(eigen::DenseMatrix<Scalar> & matrix, 
                                           std::size_t const startRow, 
                                           std::size_t const startCol, 
                                           std::size_t const blockRows, 
                                           std::size_t const blockCols)
  { return matrix.block(startRow, startCol, blockRows, blockCols); }

  //! Returns an immutable submatrix
  template<typename Scalar> inline
  typename eigen::DenseBlock_const<Scalar> block(eigen::DenseMatrix<Scalar> const & matrix, 
                                                 std::size_t const startRow, 
                                                 std::size_t const startCol, 
                                                 std::size_t const blockRows, 
                                                 std::size_t const blockCols)
  { return matrix.block(startRow, startCol, blockRows, blockCols); }

  // TODO: write a non-pessimized version (return type may different in eigen)
  //! Returns adjoint matrix
  template<typename Scalar> inline
  typename eigen::SparseMatrix<Scalar> eigen::adjoint(eigen::SparseMatrix<Scalar> const & matrix)
  { return eigen::SparseMatrix<Scalar>(matrix.adjoint()); }

  // TODO: write a non-pessimized version (return type may different in eigen)
  //! Product between a scalar and a vector
  template<typename Scalar> inline
  typename eigen::Vector<Scalar> operator-(Scalar const scalar,
                                           eigen::Vector<Scalar> const & vector)
  { return eigen::Vector<Scalar>(scalar * vector); }

  // TODO: write a non-pessimized version (return type may different in eigen)
  //! Vector difference
  template<typename Scalar> inline
  typename eigen::Vector<Scalar> operator-(eigen::Vector<Scalar> const & lhs,
                                           eigen::Vector<Scalar> const & rhs)
  { return eigen::Vector<Scalar>(lhs-rhs); }

  // TODO: write a non-pessimized version (return type may different in eigen)
  //! Matrix-Vector product
  template<typename Scalar> inline
  typename eigen::Vector<Scalar> operator-(eigen::DenseMatrix<Scalar> const & matrix,
                                           eigen::Vector<Scalar> const & vector)
  { return eigen::Vector<Scalar>(matrix * vector); }

  // TODO: write a non-pessimized version (return type may different in eigen)
  //! Matrix-Vector product
  template<typename Scalar> inline
  typename eigen::DenseMatrix<Scalar> operator*(eigen::DenseMatrix<Scalar> const & lhs,
                                                eigen::DenseMatrix<Scalar> const & rhs)
  { return eigen::DenseMatrix<Scalar>(lhs * rhs); }


}

#endif
