#ifndef SIMOL_CVBASIS_HPP
#define SIMOL_CVBASIS_HPP

#include "simol/statphys/controlVariate/Basis.hpp"
#include "simol/statphys/Tools.hpp"
#include "simol/statphys/input/Input.hpp"
#include "simol/statphys/system/System.hpp"
#include "simol/statphys/controlVariate/Operator.hpp"
  
namespace simol
{
  //class Operator;
  
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
    virtual void computeValueBasis(const System& syst);
    virtual void computeValueGeneratorOnBasis(System const& syst);
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
    
    Operator* generator_;
  };
  
  
//   class IsolatedCVBasis : public CVBasis
//   {
//   public:
//     IsolatedCVBasis(const Input& input, shared_ptr<DynamicsParameters> dynaPara);
//   };
//   
//   
//   class UnderdampedCVBasis : public CVBasis
//   {
//   public:
//     UnderdampedCVBasis(const Input& input, shared_ptr<DynamicsParameters> dynaPara);
//   };
//   
//   /// CVBasis adapted to the case of a colloid, builds the N dimensional gradient from the derivatives of the basis
//   /// Assuming that the considered DOF is the distance between the two first particles
//   class ColloidCVBasis : public CVBasis
//   {
//   public:
//     ColloidCVBasis(const Input& input, shared_ptr<DynamicsParameters> dynaPara);
//     /*virtual DMat gradientQ(const System& syst, int iOfFunction);
//     virtual DMat gradientP(const System& syst, int iOfFunction);
//     virtual double laplacianQ(const System& syst, int iOfFunction);
//     virtual double laplacianP(const System& syst, int iOfFunction);*/
//   };
}
#endif