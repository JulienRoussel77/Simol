/* 
 * File:   SchrodingerHamiltonian.hpp
 * Author: cdoucet
 *
 * Created on 28 octobre 2015, 14:19
 */

#ifndef SIMOL_SCHRODINGERHAMILTONIAN_HPP
#define	SIMOL_SCHRODINGERHAMILTONIAN_HPP

namespace simol
{
    class SchrodingerHamiltonian
    {
        public:
            SchrodingerHamiltonian(std::string const & pathToData)
            : kinetic_(MatrixMarketFile(pathToData + "kinetic_matrix.mtx")),
              overlap_(MatrixMarketFile(pathToData + "overlap_matrix.mtx")),
              potential_(MatrixMarketFile(pathToData + "potential_matrix.mtx"))
            {
                // Need to define a file format for tensors to be able to use list initializations for tensors
                twoElectrons_ = SparseTensor<double>(pathToData + "twoelectron_matrix.txt", basisDimension());
            }
            
            std::size_t basisDimension() const
            { return kinetic_.numberOfRows(); }
            
        private:
            DenseMatrix<double> kinetic_;
            DenseMatrix<double> overlap_;
            DenseMatrix<double> potential_;
            SparseTensor<double> twoElectrons_;
    };
}

#endif	/* SIMOL_SCHRODINGERHAMILTONIAN_HPP */

