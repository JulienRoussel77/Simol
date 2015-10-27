#ifndef SIMOL_SLATERDETERMINANT_HPP
#define	SIMOL_SLATERDETERMINANT_HPP

#include "../linalg/Matrix.hpp"

namespace simol
{
    class SlaterDeterminant
    {
        public:
            SlaterDeterminant(DenseMatrix<double> const & matrix)
            : matrix_(matrix)
            {}
            DenseMatrix<double> const & matrix() const
            { return matrix_; }
        private:
            DenseMatrix<double>  matrix_;
    };
    
    DenseMatrix<double>
    Smat(SlaterDeterminant const & Phi, 
         SlaterDeterminant const & Psi, 
         SparseMatrix<double> const & overlap) // mal nommé : utilisé avec H
    {
        DenseMatrix<double> U = Phi.matrix();
        DenseMatrix<double> V = Psi.matrix();

        // recupere la partie non stockee (symetrique)
        DenseMatrix<double> temp(overlap.numberOfRows(), V.number_of_columns());
        temp.wrapped_ = overlap.wrapped_.selfadjointView<Eigen::Upper>() * V.wrapped_;

        return DenseMatrix<double>(U.wrapped_.adjoint() * temp.wrapped_);
    }

    
 

}


#endif	

