#ifndef SIMOL_CVBASIS_HPP
#define SIMOL_CVBASIS_HPP

#include "simol/statphys/controlVariate/Basis.hpp"
#include "simol/statphys/Tools.hpp"
#include "simol/statphys/input/Input.hpp"
  
namespace simol
{
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