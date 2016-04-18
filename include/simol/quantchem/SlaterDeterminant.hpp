#ifndef SIMOL_SLATERDETERMINANT_HPP
#define SIMOL_SLATERDETERMINANT_HPP

#include "simol/core/linalg/DenseMatrix.hpp"

namespace simol
{
  class SlaterDeterminant
  {
    public:
      SlaterDeterminant(DenseMatrix<double> const & matrix);
      DenseMatrix<double> const & matrix() const;
      std::size_t number_of_electrons() const;
    private:
      DenseMatrix<double>  matrix_;
  };

  DenseMatrix<double>
  Smat(SlaterDeterminant const & Phi,
       SlaterDeterminant const & Psi,
       DenseMatrix<double> const & overlap); // mal nommé : utilisé avec H


}


#endif

