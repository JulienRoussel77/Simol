#ifndef SIMOL_DISCRETEHAMILTONIAN_HPP
#define	SIMOL_DISCRETEHAMILTONIAN_HPP

namespace simol
{
    class DiscreteHamiltonian
    {
        public:
           DiscreteHamiltonian(std::string pathToData, std::size_t basisDimension)
           : kinetic_(pathToData + "kinetic_matrix.txt", basisDimension),
             overlap_(pathToData + "overlap_matrix.txt", basisDimension),
             potential_(pathToData + "potential_matrix.txt", basisDimension),
             two_electrons_(pathToData + "twoelectron_matrix.txt", basisDimension)
           {}

           SparseMatrix<double> const &
           kinetic() const
           { return kinetic_; }

           SparseMatrix<double> const &
           overlap() const
           { return overlap_; }

           SparseMatrix<double> const &
           potential() const
           { return potential_; }

           SparseTensor<double> const &
           two_electrons() const
           { return two_electrons_; }

            std::size_t
            basisDimension() const
            { return kinetic_.numberOfRows(); }

        private:
            SparseMatrix<double> kinetic_;
            SparseMatrix<double> overlap_;
            SparseMatrix<double> potential_;
            SparseTensor<double> two_electrons_;
    };
}

#endif	/* SIMOL_DISCRETEHAMILTONIAN_HPP */

