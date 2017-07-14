#ifndef SIMOL_CVBASIS_HPP
#define SIMOL_CVBASIS_HPP

#include "simol/statphys/controlVariate/Basis.hpp"
#include "simol/statphys/Tools.hpp"
#include "simol/statphys/input/Input.hpp"
#include "simol/statphys/system/System.hpp"
  
namespace simol
{
  class CVBasis;
  shared_ptr<CVBasis> createCVBasis(Input const& input);
  
  /// Structure containing a basis type
  /// basisValues_ is a vector containing the values of the basis functions e_k(q,p)
  /// generatorOnBasisValues_ is a vector containing the values of the functions (L e_k)(q,p)
  class CVBasis
  {
    friend shared_ptr<CVBasis> createCVBasis(Input const& input);
  public:
    CVBasis(Input const& input);
    ~CVBasis();
    int const& totalNbOfElts() const;
    virtual DMat qVariable(const System& syst) const;
    virtual DMat pVariable(const System& syst) const;
    virtual DMat forces(const System& syst) const;
    virtual DVec basisVariables(const System& syst) const = 0;
    virtual void computeValueBasis(const System& syst);
    virtual DMat gradientQ(const System& syst, int iOfFunction);
    virtual DMat gradientP(const System& syst, int iOfFunction);
    virtual double laplacianQ(const System& syst, int iOfFunction);
    virtual double laplacianP(const System& syst, int iOfFunction);
    
    Potential* potential_;
    TensorBasis* tensorBasis_;
    int totalNbOfElts_;
    shared_ptr<DVec> cvCoeffs_;
    
    DVec basisValues_;
    DVec generatorOnBasisValues_;
  };
  
  /// CVBasis adapted to the case of a colloid, builds the N dimensional gradient from the derivatives of the basis
  /// Assuming that the considered DOF is the distance between the two first particles
  class IsolatedCVBasis : public CVBasis
  {
  public:
    IsolatedCVBasis(const Input& input);
    virtual DVec basisVariables(const System& syst) const;
  };
  
  /// CVBasis adapted to the case of a colloid, builds the N dimensional gradient from the derivatives of the basis
  /// Assuming that the considered DOF is the distance between the two first particles
  class ColloidCVBasis : public CVBasis
  {
  public:
    ColloidCVBasis(const Input& input);
    virtual DMat qVariable(const System& syst) const;
    virtual DMat pVariable(const System& syst) const;
    virtual DMat forces(const System& syst) const;
    virtual DVec basisVariables(const System& syst) const;   
    /*virtual DMat gradientQ(const System& syst, int iOfFunction);
    virtual DMat gradientP(const System& syst, int iOfFunction);
    virtual double laplacianQ(const System& syst, int iOfFunction);
    virtual double laplacianP(const System& syst, int iOfFunction);*/
  };
  
  /*class ExactColloidCVBasis : public ColloidCVBasis
  {
  public:
    double meshStep;
    double xmin;
    DVec 
  };*/
}
#endif