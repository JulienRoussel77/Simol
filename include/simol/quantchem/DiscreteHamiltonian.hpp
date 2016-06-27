#ifndef SIMOL_DISCRETEHAMILTONIAN_HPP
#define SIMOL_DISCRETEHAMILTONIAN_HPP

#include "simol/core/linalg/SparseTensor.hpp"
#include "simol/core/linalg/SymSparseMatrix.hpp"

#include <string>

namespace simol
{

  class DiscreteHamiltonian
  {
    public:
      DiscreteHamiltonian(std::string pathToData, std::size_t basisDimension);

      SymSparseMatrix<double> const & kinetic() const;
      SymSparseMatrix<double> const & overlap() const;
      SymSparseMatrix<double> const & potential() const;
      SparseTensor<double> const & two_electrons() const;

      std::size_t basisDimension() const;

    private:
      SymSparseMatrix<double> kinetic_;
      SymSparseMatrix<double> overlap_;
      SymSparseMatrix<double> potential_;
      SparseTensor<double> two_electrons_;
  };
}

#endif  /* SIMOL_DISCRETEHAMILTONIAN_HPP */

