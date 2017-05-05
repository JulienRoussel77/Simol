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
    //cout << "cvCoeffs : " << endl << reshape(*cvCoeffs0, basis_->nbOfElts(0), basis_->nbOfElts(1)) << endl;
  }
  
  CVBasis::~CVBasis()
  {
    //if (basis_) delete basis_;
  }

  const int& CVBasis::totalNbOfElts() const
  {return totalNbOfElts_;}
  
  DMat CVBasis::qVariable(System const& syst) const
  {
    DMat qVariableMat= DVec::Zero(1);
    qVariableMat(0,0) = syst(0).position(0);
    return qVariableMat;
  }
  
  DMat CVBasis::pVariable(System const& syst) const
  {
    /*DMat momentumMat = DVec::Zero(1);
    DVec r01 = syst(1).position() - syst(0).position();
    double d01 = r01.norm();
    momentumMat(0,0) = dot(syst(1).momentum() - syst(0).momentum(), r01) / d01;*/
    
    DMat qVariableMat = DVec::Zero(1);
    qVariableMat(0,0) = syst(0).momentum(0);
    return qVariableMat;
  }
  
  DMat CVBasis::forces(System const& syst) const
  {
    /*DMat forcesMat = DVec::Zero(1);
    DVec r01 = syst(1).position() - syst(0).position();
    double d01 = r01.norm();
    forcesMat(0,0) = dot(syst(1).force() - syst(0).force(), r01) / d01;*/
    
    DMat forceMat = DVec::Zero(1);
    forceMat(0,0) = syst(0).force(0);
    return forceMat;
  }
  
  ///
  ///Evaluates each function of the basis at "configuration" and updates the attribute
  void CVBasis::computeValueBasis(System const& syst)
  {
    DVec variables = basisVariables(syst);
    for (int iOfElt = 0; iOfElt < totalNbOfElts(); iOfElt++)
      basisValues_(iOfElt) = basis_->value(variables, iOfElt); 
  }
  
  DMat CVBasis::gradientQ(System const& syst, int iOfFunction)
  {
    return basis_->gradientQ(basisVariables(syst), iOfFunction);    
  }
  
  DMat CVBasis::gradientP(System const& syst, int iOfFunction)
  {
    return basis_->gradientP(basisVariables(syst), iOfFunction);    
  }
  
  double CVBasis::laplacianQ(System const& syst, int iOfFunction)
  {
    return basis_->laplacianQ(basisVariables(syst), iOfFunction);
  }
  
  double CVBasis::laplacianP(System const& syst, int iOfFunction)
  {
    return basis_->laplacianP(basisVariables(syst), iOfFunction);    
  }
  
  
  
  IsolatedCVBasis::IsolatedCVBasis(TensorBasis* basis0, shared_ptr<DVec> cvCoeffs0):
    CVBasis(basis0, cvCoeffs0)
  {}
  
  DVec IsolatedCVBasis::basisVariables(System const& syst) const
  {
    DVec variables = DVec::Zero(2);
    variables(0) = syst(0).position(0);
    variables(1) = syst(0).momentum(0);
    return variables;
  }
  


  
  ColloidCVBasis::ColloidCVBasis(TensorBasis* basis0, shared_ptr<DVec> cvCoeffs0):
    CVBasis(basis0, cvCoeffs0)
  {}
  
  /*DMat ColloidCVBasis::gradientQ(System const& syst, int iOfFunction)
  {
    return 2*basis_->gradientQ(basisVariables(syst), iOfFunction);    
  }
  
  DMat ColloidCVBasis::gradientP(System const& syst, int iOfFunction)
  {
    return 2*basis_->gradientP(basisVariables(syst), iOfFunction);    
  }
  
  double ColloidCVBasis::laplacianQ(System const& syst, int iOfFunction)
  {
    return 2*syst.dimension() * basis_->laplacianQ(basisVariables(syst), iOfFunction);
  }
  
  double ColloidCVBasis::laplacianP(System const& syst, int iOfFunction)
  {
    return 2*syst.dimension() * basis_->laplacianP(basisVariables(syst), iOfFunction);    
  }*/
  
  DVec ColloidCVBasis::basisVariables(System const& syst) const
  {
    DVec variables = DVec::Zero(2);
    DVec r01np = syst(1).position() - syst(0).position();
    DVec r01 = syst.periodicDistance(r01np);
    double d01 = r01.norm();
    DVec p01 = syst(1).momentum() - syst(0).momentum();
    variables(0) = d01;
    variables(1) = dot(p01, r01)/d01;
    return variables;
  }
  
  DMat ColloidCVBasis::qVariable(System const& syst) const
  {
    DMat distanceMat = DVec::Zero(1);
    DVec r01np = syst(1).position() - syst(0).position();
    DVec r01 = syst.periodicDistance(r01np);
    distanceMat(0,0) = r01.norm();
    return distanceMat;
  }
  
  DMat ColloidCVBasis::pVariable(System const& syst) const
  {
    DMat momentumMat = DVec::Zero(1);
    DVec r01np = syst(1).position() - syst(0).position();
    DVec r01 = syst.periodicDistance(r01np);
    double d01 = r01.norm();
    momentumMat(0,0) = dot(syst(1).momentum() - syst(0).momentum(), r01) / d01;
    return momentumMat;
  }
  
  DMat ColloidCVBasis::forces(System const& syst) const
  {
    DMat forcesMat = DVec::Zero(1);
    DVec r01np = syst(1).position() - syst(0).position();
    DVec r01 = syst.periodicDistance(r01np);
    double d01 = r01.norm();
    forcesMat(0,0) = dot(syst(1).force() - syst(0).force(), r01) / (2*d01) + (syst.dimension()-1) / d01;
    //cout << r01.adjoint() << " -> "  << (syst(1).force() - syst(0).force()).adjoint() << endl;
    //cout << d01 << " -> " << (syst(1).force() - syst(0).force()).norm() << endl;
    return forcesMat;
  }
  
}   
  