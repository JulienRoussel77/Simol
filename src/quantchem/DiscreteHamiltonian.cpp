
#include "simol/quantchem/DiscreteHamiltonian.hpp"

namespace simol
{
    DiscreteHamiltonian::DiscreteHamiltonian(std::string pathToData, std::size_t basisDimension)
    : kinetic_(pathToData + "kinetic_matrix.txt", basisDimension),
      overlap_(pathToData + "overlap_matrix.txt", basisDimension),
      potential_(pathToData + "potential_matrix.txt", basisDimension),
      two_electrons_(pathToData + "twoelectron_matrix.txt", basisDimension)
    {}

    SparseMatrix<double> const &
    DiscreteHamiltonian::kinetic() const
    { return kinetic_; }

    SparseMatrix<double> const &
    DiscreteHamiltonian::overlap() const
    { return overlap_; }

    SparseMatrix<double> const &
    DiscreteHamiltonian::potential() const
    { return potential_; }

    SparseTensor<double> const &
    DiscreteHamiltonian::two_electrons() const
    { return two_electrons_; }

     std::size_t
     DiscreteHamiltonian::basisDimension() const
     { return kinetic_.numberOfRows(); }

}
