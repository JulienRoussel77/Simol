#include "simol/statphys/controlVariate/CVBasis.hpp"
#include <simol/core/io/CommandLine.hpp>

namespace simol
{
  /*CVBasis::CVBasis(Input const& input):
    nbOfFunctions_(input.nbOfFourier()*input.nbOfHermite()),
    basis_(createBasis(input)),,
    basisValues_(input.nbOfFourier() * input.nbOfHermite()),
    generatorOnBasisValues_(input.nbOfFourier() * input.nbOfHermite())
  {}*/
  
  CVBasis::CVBasis(TensorBasis* basis0, shared_ptr<DVec> cvCoeffs0):
    totalNbOfElts_(basis0->totalNbOfElts()),
    basis_(basis0),
    cvCoeffs_(cvCoeffs0),
    basisValues_(totalNbOfElts_),
    generatorOnBasisValues_(totalNbOfElts_)
  {
    cout << "cvCoeffs : " << endl << reshape(*cvCoeffs0, basis_->nbOfElts(0), basis_->nbOfElts(1)) << endl;
  }
  
  CVBasis::~CVBasis()
  {
    //if (basis_) delete basis_;
  }

  const int& CVBasis::totalNbOfElts() const
  {return totalNbOfElts_;}
  
  ///
  ///Evaluates each function of the basis at "configuration" and updates the attribute
  void CVBasis::computeValueBasis(System const& syst)
  {
    for (int iOfElt = 0; iOfElt < totalNbOfElts(); iOfElt++)
      basisValues_(iOfElt) = basis_->value(syst, iOfElt);
    
    //ofstream tempOut("valueBasis.txt", std::ofstream::app);
    //tempOut << modulo(syst(0).position(0), 0, 2*M_PI) << " " << syst(0).force(0) << " " << dot(*cvCoeffs_, basisValues_)-1 << endl;
  }
}   
  