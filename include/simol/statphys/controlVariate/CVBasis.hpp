#ifndef SIMOL_CVBASIS_HPP
#define SIMOL_CVBASIS_HPP

#include "simol/statphys/controlVariate/Basis.hpp"
#include "simol/statphys/Tools.hpp"
#include "simol/statphys/input/Input.hpp"
#include "simol/statphys/system/System.hpp"
  
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
    virtual DMat qVariable(System const& syst) const;
    virtual DMat pVariable(System const& syst) const;
    virtual DMat forces(System const& syst) const;
    virtual DVec basisVariables(System const& syst) const = 0;
    virtual void computeValueBasis(System const& syst);
    virtual DMat gradientQ(System const& syst, int iOfFunction);
    virtual DMat gradientP(System const& syst, int iOfFunction);
    virtual double laplacianQ(System const& syst, int iOfFunction);
    virtual double laplacianP(System const& syst, int iOfFunction);
    
    int totalNbOfElts_;
    TensorBasis* basis_;
    shared_ptr<DVec> cvCoeffs_;
    
    DVec basisValues_;
    DVec generatorOnBasisValues_;
  };
  
  /// CVBasis adapted to the case of a colloid, builds the N dimensional gradient from the derivatives of the basis
  /// Assuming that the considered DOF is the distance between the two first particles
  class IsolatedCVBasis : public CVBasis
  {
  public:
    IsolatedCVBasis(TensorBasis* basis0, shared_ptr<DVec> cvCoeffs0);
    virtual DVec basisVariables(System const& syst) const;
  };
  
  /// CVBasis adapted to the case of a colloid, builds the N dimensional gradient from the derivatives of the basis
  /// Assuming that the considered DOF is the distance between the two first particles
  class ColloidCVBasis : public CVBasis
  {
  public:
    ColloidCVBasis(TensorBasis* basis0, shared_ptr<DVec> cvCoeffs0);
    virtual DMat qVariable(System const& syst) const;
    virtual DMat pVariable(System const& syst) const;
    virtual DMat forces(System const& syst) const;
    virtual DVec basisVariables(System const& syst) const;   
    /*virtual DMat gradientQ(System const& syst, int iOfFunction);
    virtual DMat gradientP(System const& syst, int iOfFunction);
    virtual double laplacianQ(System const& syst, int iOfFunction);
    virtual double laplacianP(System const& syst, int iOfFunction);*/
  };
}
#endif