#ifndef SIMOL_CVBASIS_HPP
#define SIMOL_CVBASIS_HPP

#include "simol/statphys/controlVariate/Basis.hpp"
#include "simol/statphys/Tools.hpp"
#include "simol/statphys/input/Input.hpp"
  
namespace simol
{
  /// Structure containing a basis type
  /// basisValues_ is a vector containing the values of the basis functions e_k(q,p)
  /// generatorOnBasisValues_ is a vector containing the values of the functions (L e_k)(q,p)
  class CVBasis
  {
  public:
    CVBasis();
    //CVBasis(Input const& input);
    CVBasis(TensorBasis* basis0, shared_ptr<DVec> cvCoeffs0);
    ~CVBasis();
    int const& totalNbOfElts() const;
    void computeValueBasis(System const& syst);
    
    int totalNbOfElts_;
    TensorBasis* basis_;
    shared_ptr<DVec> cvCoeffs_;
    
    DVec basisValues_;
    DVec generatorOnBasisValues_;
  };
}
#endif