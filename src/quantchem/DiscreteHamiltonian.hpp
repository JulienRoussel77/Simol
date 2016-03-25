#ifndef SIMOL_DISCRETEHAMILTONIAN_HPP
#define	SIMOL_DISCRETEHAMILTONIAN_HPP

#include "core/linalg/SparseTensor.hpp"

#include <string>

namespace simol
{

    class DiscreteHamiltonian
    {
        public:
           DiscreteHamiltonian(std::string pathToData, std::size_t basisDimension);

           SparseMatrix<double> const & kinetic() const;
           SparseMatrix<double> const & overlap() const;
           SparseMatrix<double> const & potential() const;
           SparseTensor<double> const & two_electrons() const;

           std::size_t basisDimension() const;

        private:
            SparseMatrix<double> kinetic_;
            SparseMatrix<double> overlap_;
            SparseMatrix<double> potential_;
            SparseTensor<double> two_electrons_;
    };
}

#endif	/* SIMOL_DISCRETEHAMILTONIAN_HPP */

