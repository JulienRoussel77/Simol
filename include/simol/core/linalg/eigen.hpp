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

  template<typename Scalar>
  struct eigen
  {
    using DenseMatrix      = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using DenseAdjoint     = typename DenseMatrix::AdjointReturnType;
    using DenseColumn      = typename DenseMatrix::ColXpr;
    using DenseBlock       = Eigen::Block<DenseMatrix>;
    using DenseBlock_const = Eigen::Block<DenseMatrix const>;
    using Vector           = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using VectorMap        = Eigen::Map<Vector>;
    using SparseMatrix     = Eigen::SparseMatrix<Scalar>;
    using SVD              = Eigen::JacobiSVD<DenseMatrix>;

    static void fill(DenseMatrix const & matrix, Scalar scalar);

    static std::size_t size(Vector const & vector);
    static std::size_t min_index(Vector const & vector);
    static std::size_t number_of_rows(DenseMatrix const & matrix);
    static std::size_t number_of_rows(SparseMatrix const & matrix);
    static std::size_t number_of_columns(DenseMatrix const & matrix);
    static std::size_t number_of_columns(SparseMatrix const & matrix);

    static Scalar min(Vector const & vector);
    static Scalar max(Vector const & vector);
    static Scalar norm(Vector const & vector);
    static Scalar inner_product(Vector const & lhs, Vector const & rhs);
    static Scalar determinant(DenseMatrix const & matrix);
    static Scalar trace(DenseMatrix const & matrix);

    static Vector sort(Vector const & vector);
    static Vector subvector(Vector const & vector, std::size_t start, std::size_t length);
    static Vector Zero(std::size_t const length);

    static DenseMatrix Zero(std::size_t const numberOfRows, std::size_t const numberOfColumns);
    static DenseMatrix Identity(std::size_t const dimension);
    static DenseMatrix inverse(DenseMatrix const & matrix);

    static SparseMatrix adjoint(SparseMatrix const & matrix);
  };

  //! Returns the size of a vector
  template<typename Scalar> inline
  std::size_t eigen<Scalar>::size(eigen<Scalar>::Vector const & vector)
  { return vector.size(); }

  //! Returns the index of the minimum coefficient of a vector
  template<typename Scalar> inline
  std::size_t eigen<Scalar>::min_index(eigen<Scalar>::Vector const & vector)
  {
    std::size_t index;
    vector.minCoeff(&index);
    return index;
  }

  //! Returns the number of rows of a dense matrix
  template<typename Scalar> inline
  std::size_t eigen<Scalar>::number_of_rows(eigen<Scalar>::DenseMatrix const & matrix)
  { return matrix.rows(); }

  //! Returns the number of rows of a sparse matrix
  template<typename Scalar> inline
  std::size_t eigen<Scalar>::number_of_rows(eigen<Scalar>::SparseMatrix const & matrix)
  { return matrix.rows(); }

  //! Returns the number of columns of a dense matrix
  template<typename Scalar> inline
  std::size_t eigen<Scalar>::number_of_columns(eigen<Scalar>::DenseMatrix const & matrix)
  { return matrix.cols(); }

  //! Returns the number of columns of a sparse matrix
  template<typename Scalar> inline
  std::size_t eigen<Scalar>::number_of_columns(eigen<Scalar>::SparseMatrix const & matrix)
  { return matrix.cols(); }

  //! Returns the minimum coefficient of a vector
  template<typename Scalar> inline
  Scalar eigen<Scalar>::min(eigen<Scalar>::Vector const & vector)
  { return vector.minCoeff(); }

  //! Returns the maximum coefficient of a vector
  template<typename Scalar> inline
  Scalar eigen<Scalar>::max(eigen<Scalar>::Vector const & vector)
  { return vector.maxCoeff(); }

  //! Returns the norm of a vector
  template<typename Scalar> inline
  Scalar eigen<Scalar>::norm(eigen<Scalar>::Vector const & vector)
  { return vector.norm(); }

  //! Return the trace of a dense matrix
  template<typename Scalar> inline
  Scalar eigen<Scalar>::trace(eigen<Scalar>::DenseMatrix const & matrix)
  { return matrix.trace(); }

  //! Returns the determinant of a dense matrix
  template<typename Scalar> inline
  Scalar eigen<Scalar>::determinant(eigen<Scalar>::DenseMatrix const & matrix)
  { return matrix.determinant(); }


  //! Returns a sorted copy of a vector
  template<typename Scalar> inline
  typename eigen<Scalar>::Vector eigen<Scalar>::sort(eigen<Scalar>::Vector const & vector)
  {
    typename eigen<Scalar>::Vector to_be_sorted = vector;
    std::sort( to_be_sorted.data(), to_be_sorted.data() + to_be_sorted.size() );
    return to_be_sorted;
  }

  //! Returns the null vector
  template<typename Scalar> inline
  typename eigen<Scalar>::Vector eigen<Scalar>::Zero(std::size_t const length)
  { return eigen<Scalar>::Vector(eigen<Scalar>::Vector::Zero(length)); }

  // TODO: could be optimized by calling a suitable function of Eigen
  //! Returns a subvector of a vector
  template<typename Scalar>
  typename eigen<Scalar>::Vector eigen<Scalar>::subvector(eigen<Scalar>::Vector const & vector,
      std::size_t start,
      std::size_t length)
  {
    eigen<Scalar>::Vector subvec(length);
    for (std::size_t index = 0; index < length; ++index)
      subvec(index) = vector(index + start);
    return subvec;
  }

  //! Returns the inner product between two vectors
  template<typename Scalar> inline
  Scalar eigen<Scalar>::inner_product(eigen<Scalar>::Vector const & lhs,
                                      eigen<Scalar>::Vector const & rhs)
  { return lhs.dot(rhs); }

  //! Returns a null dense matrix
  template<typename Scalar> inline
  typename eigen<Scalar>::DenseMatrix eigen<Scalar>::Zero(std::size_t const numberOfRows, size_t const numberOfColumns)
  { return eigen<Scalar>::DenseMatrix(eigen<Scalar>::DenseMatrix::Zero(numberOfRows, numberOfColumns)); }

  //! Returns an identity dense matrix
  template<typename Scalar> inline
  typename eigen<Scalar>::DenseMatrix eigen<Scalar>::Identity(std::size_t const dimension)
  { return eigen<Scalar>::DenseMatrix(eigen<Scalar>::DenseMatrix::Identity(dimension)); }

  //! Returns the inverse of a dense matrix
  template<typename Scalar> inline
  typename eigen<Scalar>::DenseMatrix eigen<Scalar>::inverse(eigen<Scalar>::DenseMatrix const & matrix)
  { return eigen<Scalar>::DenseMatrix(matrix.inverse()); }

  // TODO: write a non-pessimized version (return type may different in eigen)
  //! Returns adjoint matrix
  template<typename Scalar> inline
  typename eigen<Scalar>::SparseMatrix eigen<Scalar>::adjoint(eigen<Scalar>::SparseMatrix const & matrix)
  { return eigen<Scalar>::SparseMatrix(matrix.adjoint()); }


}

#endif
