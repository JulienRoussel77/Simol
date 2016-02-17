#ifndef SIMOL_SLATERDETERMINANT_HPP
#define	SIMOL_SLATERDETERMINANT_HPP

#include "core/linalg/DenseMatrix.hpp"

namespace simol
{
    class SlaterDeterminant
    {
        public:
            SlaterDeterminant(DenseMatrix<double> const & matrix);
            DenseMatrix<double> const & matrix() const;
        private:
            DenseMatrix<double>  matrix_;
    };

    SlaterDeterminant::SlaterDeterminant(DenseMatrix<double> const & matrix)
    : matrix_(matrix)
    {}

    DenseMatrix<double> const &
    SlaterDeterminant::matrix() const
    { return matrix_; }

    DenseMatrix<double>
    Smat(SlaterDeterminant const & Phi,
         SlaterDeterminant const & Psi,
         DenseMatrix<double> const & overlap) // mal nommé : utilisé avec H
    {
        DenseMatrix<double> U = Phi.matrix();
        DenseMatrix<double> V = Psi.matrix();

        // recupere la partie non stockee (symetrique)
        DenseMatrix<double> temp(overlap.numberOfRows(), V.numberOfColumns());
        temp.wrapped_ = overlap.wrapped_.selfadjointView<Eigen::Upper>() * V.wrapped_;

        return DenseMatrix<double>(U.wrapped_.adjoint() * temp.wrapped_);
    }




}


#endif

