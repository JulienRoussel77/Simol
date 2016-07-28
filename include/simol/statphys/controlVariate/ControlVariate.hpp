#ifndef SIMOL_CONTROLVARIATE_HPP
#define SIMOL_CONTROLVARIATE_HPP

#include "simol/statphys/output/Observable.hpp"
#include "simol/statphys/Tools.hpp"
#include "simol/statphys/output/Statistics.hpp"
#include "simol/statphys/potential/Potential.hpp"
#include "simol/statphys/system/Particle.hpp"
#include "simol/statphys/controlVariate/Basis.hpp"
#include "simol/statphys/controlVariate/Galerkin.hpp"


namespace simol
{

  //class ControlVariate;
  Observable* createControlVariate(const Input& input, const string& outPath, Potential& potential, Galerkin* galerkin);

  class ControlVariate : public Observable
  {
  friend Observable* createControlVariate(const Input& input, const string& outPath, Potential& potential, Galerkin* galerkin);
  //protected:
  public:
      int dimension_;
      //int decorrelationNbOfSteps_;
      //double timeStep_;
      int nbOfFunctions_;
      int nbOfFunctionPairs_;
      //int printPeriodNbOfSteps_;
      //int nbOfAutocoPts_;

      //AutocorrelationStats statsObservable_;
      AutocorrelationStats autocoStatsBetter_;
      
      Statistics statsGeneratorOnBasis_;

      Statistics statsB1_;
      CorrelationStats statsB2_;
      Statistics statsD_;
      Vector<double> lastA_;

      Potential* potential_;
      
      Vector<double> valueBasis_; // contains the image of the basis by the generator (Phi_i(x))_i
      Vector<double> valueGeneratorOnBasis_; // contains the image of the basis by the generator (L Phi_i(x))_i
    public:
      ControlVariate(const Input& input, const string& outPath, Potential& potential, int nbOfFunctions);

      // ACCESSEURS
      virtual bool isNone() const;
      double potential(Vector<double> const& position) const;
      Vector<double> potentialDerivative(Vector<double> const& position) const;
      double potentialLaplacian(Vector<double> const& position) const;
      //const double& timeStep() const;
      //const int& decorrelationNbOfSteps() const;
      //double decorrelationTime() const;
      virtual int nbOfFunctions() const;
      virtual int nbOfFunctionPairs() const;
      //bool doOutput(long int iOfStep) const;
      //virtual int const& nbOfAutocoPts() const;
      virtual int nbOfFourier() const;
      virtual int nbOfHermite() const;

      virtual Vector<double> lastB1() const;
      virtual Vector<double> meanB1() const;
      virtual Vector<double> meanB() const;
      virtual DenseMatrix<double> lastD() const;
      virtual DenseMatrix<double> meanD() const;

      /*virtual double lastObservable() const;
      virtual double meanObservable() const;
      virtual double varObservable() const;
      virtual double stdDevObservable() const;
      virtual double varOfVarObservable() const;
      virtual double stdDevOfVarObservable() const;
      //virtual double stdErrorOfVarObservable() const;*/
      
      virtual double lastValueBetter() const;
      virtual double meanBetter() const;
      virtual double varBetter() const;
      virtual double stdDevBetter() const;
      virtual double varOfVarBetter() const;
 
      virtual Vector<double> lastGeneratorOnBasis() const;
      virtual Vector<double> meanGeneratorOnBasis() const;

      /*virtual double correlationAtSpan(int iOfSpan) const;
      virtual double unbiasedCorrelationAtSpan(int iOfSpan) const;
      virtual double stdDevCorrelationAtSpan(int iOfSpan) const;
      virtual double varCorrelationAtSpan(int iOfSpan) const;*/
      //virtual double stdErrorCorrelationAtSpan(int iOfSpan) const;
      
      virtual double autocorrelationB2(int iOfSpan, int iOfFunction = 0) const;

      virtual Vector<double> lastA() const;
      virtual double lastA(int iOfFunction = 0) const;

      virtual Vector<double> correlationB2() const;
      virtual double correlationB2(int iOfFunction) const;

      // APPEND

      void appendToObservable(double value, long int iOfStep);
      void appendToB1(double value);
      void appendToB2(double value, long int iOfStep);
      void appendToD();
      void appendToBetter(double value, long int iOfStep);

      virtual void append(double value, long int iOfStep);

      // FUNCTION CARACTERIZATION

      virtual double basisFunction(vector<Particle> const& configuration, int iOfFunction = 0) const = 0;


      virtual double laplacianQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const = 0;
      virtual Vector<double> gradientQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const = 0;
      virtual double laplacianP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const = 0;
      virtual Vector<double> gradientP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const = 0;
      
      void updateHamiltonian(vector<Particle> const& configuration);
      void updateOverdamped(vector<Particle> const& configuration, double beta);
      void updateLangevin(vector<Particle> const& configuration, double beta, double gamma);
      void updateBoundaryLangevin(vector<Particle> const& configuration, double betaLeft, double betaRight, double gamma);

      void computeValueBasis(vector<Particle> const& configuration);
      void computeGeneratorHamiltonian(vector<Particle> const& configuration);
      void computeGeneratorOverdamped(vector<Particle> const& configuration, double beta);
      void computeGeneratorLangevin(vector<Particle> const& configuration, double beta, double gamma);
      void computeGeneratorBoundaryLangevin(vector<Particle> const& configuration, double betaLeft, double betaRight, double gamma);



      virtual void display(std::ofstream& out, double time) const;

      // TEMP

      virtual void displayMap(ofstream& /*out*/) const {};
      virtual void displayGradQMap(ofstream& /*out*/) const {};
      virtual void displayGradPMap(ofstream& /*out*/) const {};
  };

  /*class NoControlVariate : public ControlVariate
  {
    public:
      NoControlVariate(const Input& input, Potential& potential);
      virtual int nbOfFunctions() const;
      virtual int nbOfFunctionPairs() const;
      bool isNone() const;
      Vector<double> lastGeneratorOnBasis() const;
      double basisFunction(vector<Particle> const& configuration, int iOfFunction = 0) const;
      virtual double laplacianQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
      virtual Vector<double> gradientQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
      virtual double laplacianP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
      virtual Vector<double> gradientP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
      void update(double value, Vector<double>& generatorOnBasisFunction, vector<Particle> const& configuration, long int iOfStep);
      virtual void postTreat(std::ofstream& out, double timeStep);
  };*/


  class SinusControlVariate : public ControlVariate
  {
    public:
      SinusControlVariate(const Input& input, const string& outPath, Potential& potential);
      double basisFunction(vector<Particle> const& configuration, int iOfFunction = 0) const;
      virtual double laplacianQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
      virtual Vector<double> gradientQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
      virtual double laplacianP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
      virtual Vector<double> gradientP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
  };

  class CosControlVariate : public ControlVariate
  {
    public:
      CosControlVariate(const Input& input, const string& outPath, Potential& potential);
      double basisFunction(vector<Particle> const& configuration, int iOfFunction = 0) const;
      virtual double laplacianQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
      virtual Vector<double> gradientQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
      virtual double laplacianP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
      virtual Vector<double> gradientP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
  };

  class SinExpControlVariate : public ControlVariate
  {
    public:
      SinExpControlVariate(const Input& input, const string& outPath, Potential& potential);
      double basisFunction(vector<Particle> const& configuration, int iOfFunction = 0) const;
      virtual double laplacianQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
      virtual Vector<double> gradientQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
      virtual double laplacianP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
      virtual Vector<double> gradientP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
  };

  class CosExpControlVariate : public ControlVariate
  {
    public:
      CosExpControlVariate(const Input& input, const string& outPath, Potential& potential);
      double basisFunction(vector<Particle> const& configuration, int iOfFunction = 0) const;
      virtual double laplacianQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
      virtual Vector<double> gradientQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
      virtual double laplacianP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
      virtual Vector<double> gradientP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
  };

  class LangevinControlVariate : public ControlVariate
  {
    public:
      LangevinControlVariate(const Input& input, const string& outPath, Potential& potential);
      double basisFunction(vector<Particle> const& configuration, int iOfFunction = 0) const;
      virtual double laplacianQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
      virtual Vector<double> gradientQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
      virtual double laplacianP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
      virtual Vector<double> gradientP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
  };

  class SumEnergyControlVariate : public ControlVariate
  {
      double i0_;
    public:
      SumEnergyControlVariate(const Input& input, const string& outPath, Potential& potential);
      double basisFunction(vector<Particle> const& configuration, int iOfFunction = 0) const;
      virtual double laplacianQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
      virtual Vector<double> gradientQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
      virtual double laplacianP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
      virtual Vector<double> gradientP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
  };

  class EnergyControlVariate : public ControlVariate
  {
      double i0_;
    public:
      EnergyControlVariate(const Input& input, const string& outPath, Potential& potential);
      double basisFunction(vector<Particle> const& configuration, int iOfFunction = 0) const;
      virtual double laplacianQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
      virtual Vector<double> gradientQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
      virtual double laplacianP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
      virtual Vector<double> gradientP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
  };

  class LocalControlVariate : public ControlVariate
  {
    public:
      LocalControlVariate(const Input& input, const string& outPath, Potential& potential);
      double basisFunction(vector<Particle> const& configuration, int iOfFunction = 0) const;
      virtual double laplacianQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
      virtual Vector<double> gradientQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
      virtual double laplacianP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
      virtual Vector<double> gradientP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
  };

  class KineticControlVariate : public ControlVariate
  {
    public:
      KineticControlVariate(const Input& input, const string& outPath, Potential& potential);
      double basisFunction(vector<Particle> const& configuration, int iOfFunction = 0) const;
      virtual double laplacianQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
      virtual Vector<double> gradientQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
      virtual double laplacianP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
      virtual Vector<double> gradientP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
  };



  class TwoControlVariate : public ControlVariate
  {
    public:
      TwoControlVariate(const Input& input, const string& outPath, Potential& potential);
      double basisFunction(vector<Particle> const& configuration, int iOfFunction = 0) const;
      virtual double laplacianQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
      virtual Vector<double> gradientQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
      virtual double laplacianP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
      virtual Vector<double> gradientP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
  };

  class BasisControlVariate : public ControlVariate
  {
    protected:
      SMat coeffsVec_;
      TensorBasis* basis_;
    public:
      BasisControlVariate(const Input& input, const string& outPath, Potential& potential, Galerkin* galerkin);
      double basisFunction(vector<Particle> const& configuration, int iOfFunction = 0) const;
      virtual double laplacianQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
      virtual Vector<double> gradientQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
      virtual double laplacianP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
      virtual Vector<double> gradientP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
  };

  class ExpFourierHermiteControlVariate : public BasisControlVariate
  {
      int nbQ_, nbP_;
      double pMax_, deltaQ_, deltaP_;
    public:
      ExpFourierHermiteControlVariate(const Input& input, const string& outPath, Potential& potential, Galerkin* galerkin);
      int nbOfFourier() const;
      int nbOfHermite() const;
      virtual void displayMap(ofstream& out) const;
      virtual void displayGradQMap(ofstream& out) const;
      virtual void displayGradPMap(ofstream& out) const;
  };

}

#endif
