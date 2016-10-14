#ifndef SIMOL_MATRIXFREE_HPP
#define SIMOL_MATRIXFREE_HPP

#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>
#include "simol/statphys/Tools.hpp"
class MatrixFunction;
using Eigen::SparseMatrix;
namespace Eigen {
namespace internal {
  // MatrixFunction looks-like a SparseMatrix, so let's inherits its traits:
  template<>
  struct traits<MatrixFunction> :  public internal::traits<SMat >
  {};
}
}
// Example of a matrix-free wrapper from a user type to Eigen's compatible type
// For the sake of simplicity, this example simply wrap a Eigen::SparseMatrix.
class MatrixFunction : public Eigen::EigenBase<MatrixFunction> {
public:
  // Required typedefs, constants, and method:
  typedef double Scalar;
  typedef double RealScalar;
  typedef int StorageIndex;
  enum {
    ColsAtCompileTime = Eigen::Dynamic,
    MaxColsAtCompileTime = Eigen::Dynamic,
    IsRowMajor = false
  };
  Index rows() const { return A_->rows(); }
  Index cols() const { return A_->cols(); }
  template<typename Rhs>
  Eigen::Product<MatrixFunction,Rhs,Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs>& x) const {
    return Eigen::Product<MatrixFunction,Rhs,Eigen::AliasFreeProduct>(*this, x.derived());
  }
  // Custom API:
  MatrixFunction() : A_(0) {}
  void attachMyMatrix(const SMat &mat, DVec const& vec) {
    A_ = &mat;
    u_ = &vec;
  }
  const SMat my_matrix() const { return *A_; }
  const DVec my_vector() const { return *u_; }
private:
  const SMat *A_;
  const DVec* u_;
};
// Implementation of MatrixFunction * Eigen::DenseVector though a specialization of internal::generic_product_impl:
namespace Eigen {
namespace internal {
  template<typename Rhs>
  struct generic_product_impl<MatrixFunction, Rhs, SparseShape, DenseShape, GemvProduct> // GEMV stands for matrix-vector
  : generic_product_impl_base<MatrixFunction,Rhs,generic_product_impl<MatrixFunction,Rhs> >
  {
    typedef typename Product<MatrixFunction,Rhs>::Scalar Scalar;
    template<typename Dest>
    static void scaleAndAddTo(Dest& dst, const MatrixFunction& lhs, const Rhs& rhs, const Scalar& alpha)
    {
      // This method should implement "dst += alpha * lhs * rhs" inplace,
      // however, for iterative solvers, alpha is always equal to 1, so let's not bother about it.
      assert(alpha==Scalar(1) && "scaling is not implemented");
      dst.noalias() += lhs.my_matrix() * rhs + dot(lhs.my_vector(), rhs) * rhs;
      
      // Here we could simply call dst.noalias() += lhs.my_matrix() * rhs,
      // but let's do something fancier (and less efficient):
      for(Index i=0; i<lhs.cols(); ++i)
        dst += rhs(i) * lhs.my_matrix().col(i);
    }
  };
}
}

#endif