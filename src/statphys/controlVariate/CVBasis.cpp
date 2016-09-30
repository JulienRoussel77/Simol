#include "simol/statphys/controlVariate/CVBasis.hpp"
#include <simol/core/io/CommandLine.hpp>

namespace simol
{
  CVBasis::CVBasis(Input const& input):
    nbOfFunctions_(input.nbOfFourier()*input.nbOfHermite()),
    basisValues_(input.nbOfFourier(), input.nbOfHermite()),
    generatorOnBasisValues_(input.nbOfFourier(), input.nbOfHermite())
  {}

  const int& CVBasis::nbOfFunctions() const
  {return nbOfFunctions_;}
  
  ///
  ///Evaluates each function of the basis at "configuration" and updates the attribute
  void CVBasis::computeValueBasis(vector<Particle*> const& configuration)
  {
    for (int iOfFunction = 0; iOfFunction < nbOfFunctions_; iOfFunction++)
      basisValues_(iOfFunction) = basis_->value(configuration, iOfFunction);
  }
}   
  