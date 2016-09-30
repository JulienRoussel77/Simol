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
    CVBasis(Input const& input);
    int const& nbOfFunctions() const;
    void computeValueBasis(vector<Particle*> const& configuration);
    
    int nbOfFunctions_;
    TensorBasis* basis_;
    DVec basisValues_;
    DVec generatorOnBasisValues_;
  };
}
#endif