#ifndef SIMOL_EIGENDECOMPOSITION_HPP
#define SIMOL_EIGENDECOMPOSITION_HPP

#include "eigen.hpp"


namespace simol
{
    template<typename Scalar>
    class EigenDecomposition
    {
            typedef typename eigen<Scalar>::DenseMatrixType WrappedMatrix;
            typedef Eigen::GeneralizedSelfAdjointEigenSolver<WrappedMatrix> EigenSolver;

        public:
            EigenDecomposition(DenseMatrix<Scalar, eigen> const & leftMatrix,
                               DenseMatrix<Scalar, eigen> const & rightMatrix)
            : wrapped_(leftMatrix.wrapped_,
                       rightMatrix.wrapped_,
                       Eigen::ComputeEigenvectors|Eigen::Ax_lBx)
            {}

            Vector<double> eigenvalues()
            { return wrapped_.eigenvalues(); }
            DenseMatrix<double> eigenvectors()
            { return wrapped_.eigenvectors(); }
        private:
            EigenSolver wrapped_;
    };

}


#endif
