#include "simol/statphys/controlVariate/CVBasis.hpp"
#include <simol/core/io/CommandLine.hpp>

namespace simol
{
  CVBasis::CVBasis(Input const& input):
    nbOfFunctions_(input.nbOfFourier()*input.nbOfHermite()),
    basis_(nullptr),
    basisValues_(input.nbOfFourier() * input.nbOfHermite()),
    generatorOnBasisValues_(input.nbOfFourier() * input.nbOfHermite())
  {}
  
  CVBasis::~CVBasis()
  {
    cout << "Destruction CVBasis" << endl;
    if (basis_) delete basis_;
    cout << "OK" << endl;
  }

  const int& CVBasis::nbOfFunctions() const
  {return nbOfFunctions_;}
  
  ///
  ///Evaluates each function of the basis at "configuration" and updates the attribute
  void CVBasis::computeValueBasis(System const& syst)
  {
    for (int iOfFunction = 0; iOfFunction < nbOfFunctions_; iOfFunction++)
      basisValues_(iOfFunction) = basis_->value(syst, iOfFunction);
  }
}   
  