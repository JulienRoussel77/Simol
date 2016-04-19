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

  template<class Scalar>
  struct eigen
  {
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
    typedef Eigen::SparseMatrix<Scalar> SparseMatrixType;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> DenseMatrixType;
    typedef Eigen::Map<Vector> VectorMap;
    typedef Eigen::JacobiSVD<DenseMatrixType> SVDType;
    typedef typename DenseMatrixType::AdjointReturnType AdjointReturnType;
    typedef typename DenseMatrixType::ColXpr DenseMatrixColumn;
    typedef Eigen::Block<DenseMatrixType> DenseBlock;
    typedef Eigen::Block<DenseMatrixType const> DenseBlock_const;

    static std::size_t size(Vector const & vector);
    static std::size_t min_index(Vector const & vector);
    static Scalar min(Vector const & vector);
    static Scalar max(Vector const & vector);
    static Scalar norm(Vector const & vector);
    static Vector sort(Vector const & vector);
  };

  template<typename Scalar> inline
  std::size_t eigen<Scalar>::size(eigen<Scalar>::Vector const & vector)
  { return vector.size(); }

  template<typename Scalar> inline
  std::size_t eigen<Scalar>::min_index(eigen<Scalar>::Vector const & vector)
  {
    std::size_t index;
    vector.minCoeff(&index);
    return index;
  }

  template<typename Scalar> inline
  Scalar eigen<Scalar>::min(eigen<Scalar>::Vector const & vector)
  { return vector.minCoeff(); }

  template<typename Scalar> inline
  Scalar eigen<Scalar>::max(eigen<Scalar>::Vector const & vector)
  { return vector.maxCoeff(); }

  template<typename Scalar> inline
  Scalar eigen<Scalar>::norm(eigen<Scalar>::Vector const & vector)
  { return vector.norm(); }

  template<typename Scalar> inline
  typename eigen<Scalar>::Vector eigen<Scalar>::sort(eigen<Scalar>::Vector const & vector)
  {
    typename eigen<Scalar>::Vector to_be_sorted = vector;
    std::sort( to_be_sorted.data(), to_be_sorted.data() + to_be_sorted.size() );
    return to_be_sorted;
  }

}

#endif
