#include "simol/quantchem/SlaterDeterminant.hpp"

namespace simol
{
  SlaterDeterminant::SlaterDeterminant(DenseMatrix<double> const & matrix)
    : matrix_(matrix)
  {}

  std::size_t SlaterDeterminant::number_of_electrons() const
  { return matrix_.numberOfColumns(); }

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
    temp = overlap * V;

    return DenseMatrix<double>(U.adjoint() * temp);
  }

}
