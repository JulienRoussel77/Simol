#ifndef SIMOL_CONTROLVARIATE_HPP
#define SIMOL_CONTROLVARIATE_HPP

#include "simol/statphys/output/Observable.hpp"
#include "simol/statphys/Tools.hpp"
#include "simol/statphys/output/Statistics.hpp"
#include "simol/statphys/potential/Potential.hpp"
#include "simol/statphys/system/Particle.hpp"
#include "simol/statphys/controlVariate/Basis.hpp"
#include "simol/statphys/controlVariate/Galerkin.hpp"
#include "simol/statphys/controlVariate/CVBasis.hpp"


namespace simol
{

  //class ControlVariate;
  Observable* createControlVariate(const Input& input, int idObs, shared_ptr<CVBasis> cvBasis);

  class ControlVariate : public Observable
  {
  friend Observable* createControlVariate(const Input& input, int idObs, shared_ptr<CVBasis> cvBasis0);
  //protected:
  public:
      int dimension_;

      int nbOfFunctions_;
      int nbOfFunctionPairs_;
      
      bool doEstimateCvCoeffs_;
      
      AutocorrelationStats autocoStatsBetter_;
      
      Statistics statsGeneratorOnBasis_;

      Statistics statsB1_;
      CorrelationStats statsB2_;
      Statistics statsD_;
      DVec lastA_;

      ///
      ///cvBasis contains the basis, the cvCoeffs if they are known, and temporary basisValues
      ///it can be shared with other ControlVariates
      shared_ptr<CVBasis> cvBasis_;
    public:
      ControlVariate(const Input& input, int idObs, shared_ptr<CVBasis> cvBasis0, int nbOfFunctions);

      // ACCESSEURS
      virtual bool isNone() const;
      bool doEstimateCvCoeffs() const;
      double potential(DVec const& position) const;
      DVec potentialDerivative(DVec const& position) const;
      double potentialLaplacian(DVec const& position) const;
      virtual int nbOfFunctions() const;
      virtual int nbOfFunctionPairs() const;
      virtual int nbOfFourier() const;
      virtual int nbOfHermite() const;
      CVBasis const& cvBasis() const;
      CVBasis& cvBasis();

      virtual DVec lastB1() const;
      virtual DVec meanB1() const;
      virtual DVec meanB() const;
      virtual DMat lastD() const;
      virtual DMat meanD() const;

      
      virtual double lastValueBetter() const;
      virtual double meanBetter() const;
      virtual double varBetter() const;
      virtual double stdDevBetter() const;
      virtual double varOfVarBetter() const;
 
      virtual DVec lastGeneratorOnBasis() const;
      virtual DVec meanGeneratorOnBasis() const;
      
      virtual double autocorrelationB2(int iOfSpan, int iOfFunction = 0) const;

      virtual DVec lastA() const;
      virtual double lastA(int iOfFunction = 0) const;

      virtual DVec correlationB2() const;
      virtual double correlationB2(int iOfFunction) const;
      
      virtual double correlationBetterAtSpan(int iOfSpan) const;
      virtual double unbiasedCorrelationBetterAtSpan(int iOfSpan) const;
      virtual double varCorrelationBetterAtSpan(int iOfSpan) const;
      virtual double stdDevCorrelationBetterAtSpan(int iOfSpan) const;
      
      const double& basisValue(int iOfFunction) const;
      const DVec& basisValues() const;
      const double& generatorOnBasisValue(int iOfFunction) const;
      const DVec& generatorOnBasisValues() const;

      // APPEND

      void appendToObservable(double value, long int iOfStep);
      void appendToB1(double value);
      void appendToB2(double value, long int iOfStep);
      void appendToD();
      void appendToBetter(double value, long int iOfStep);

      virtual void append(double value, long int iOfStep);

      // FUNCTION CARACTERIZATION

      virtual double value(System const& syst) const = 0;


      /*virtual double laplacianQ(System const& syst, int iOfParticle = 0, int iOfFunction = 0) const = 0;
      virtual DVec gradientQ(System const& syst, int iOfParticle = 0, int iOfFunction = 0) const = 0;
      virtual double laplacianP(System const& syst, int iOfParticle = 0, int iOfFunction = 0) const = 0;
      virtual DVec gradientP(System const& syst, int iOfParticle = 0, int iOfFunction = 0) const = 0;*/
      

      virtual void display(long int iOfStep);
      void displayFinalValues(ofstream& out);
      void displayCorrelations(long int iOfStep);
  
      // TEMP

      /*virtual void displayMap(ofstream&) const {};
      virtual void displayGradQMap(ofstream&) const {};
      virtual void displayGradPMap(ofstream&) const {};*/
  };


  class SinusControlVariate : public ControlVariate
  {
    public:
      SinusControlVariate(const Input& input, int idObs, shared_ptr<CVBasis> cvBasis0);
      double value(System const& syst) const;
      /*virtual double laplacianQ(System const& syst, int iOfParticle = 0, int iOfFunction = 0) const;
      virtual DVec gradientQ(System const& syst, int iOfParticle = 0, int iOfFunction = 0) const;
      virtual double laplacianP(System const& syst, int iOfParticle = 0, int iOfFunction = 0) const;
      virtual DVec gradientP(System const& syst, int iOfParticle = 0, int iOfFunction = 0) const;*/
  };

  class BasisControlVariate : public ControlVariate
  {
    protected:

      //TensorBasis* basis_;
    public:
      BasisControlVariate(const Input& input, int idObs, shared_ptr<CVBasis> cvBasis0);
      DVec const& cvCoeffs() const;
      DVec& cvCoeffs();
      double const& cvCoeffs(int i) const;
      double& cvCoeffs(int i);
      
      double value(System const& syst) const;
      virtual void display(long int iOfStep);
      /*virtual double laplacianQ(System const& syst, int iOfParticle = 0, int iOfFunction = 0) const;
      virtual DVec gradientQ(System const& syst, int iOfParticle = 0, int iOfFunction = 0) const;
      virtual double laplacianP(System const& syst, int iOfParticle = 0, int iOfFunction = 0) const;
      virtual DVec gradientP(System const& syst, int iOfParticle = 0, int iOfFunction = 0) const;*/
  };

  class ExpFourierHermiteControlVariate : public BasisControlVariate
  {
      int nbQ_, nbP_;
      double pMax_, deltaQ_, deltaP_;
    public:
      ExpFourierHermiteControlVariate(const Input& input, int idObs, shared_ptr<CVBasis> cvBasis0);
      int nbOfFourier() const;
      int nbOfHermite() const;
      /*virtual void displayMap(ofstream& out) const;
      virtual void displayGradQMap(ofstream& out) const;
      virtual void displayGradPMap(ofstream& out) const;*/
  };
  
  class HermiteHermiteControlVariate : public BasisControlVariate
  {
    public:
      HermiteHermiteControlVariate(const Input& input, int idObs, shared_ptr<CVBasis> cvBasis0);
  };
  
  class ExpHermiteHermiteControlVariate : public BasisControlVariate
  {
    public:
      ExpHermiteHermiteControlVariate(const Input& input, int idObs, shared_ptr<CVBasis> cvBasis0);
  };
}

#endif
