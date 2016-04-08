#ifndef SIMOL_CONTROLVARIATE_HPP
#define SIMOL_CONTROLVARIATE_HPP

#include "Tools.hpp"
#include "Statistics.hpp"
#include "Potential.hpp"
#include "Particle.hpp"
#include "Basis.hpp"
#include "Galerkin.hpp"


namespace simol
{

  class ControlVariate;
  ControlVariate* createControlVariate(const simol::Input& input, simol::Potential& potential, Galerkin* galerkin);

  class ControlVariate
  {
  friend ControlVariate* createControlVariate(Input const& input, Potential& potential);
  protected:
    int dimension_;
    int decorrelationNbOfIterations_;
    double decorrelationTime_;
    int nbOfFunctions_;
    int nbOfFunctionPairs_;
		int periodNbOfIterations_;
		int nbOfAutocoPts_;

    AutocorrelationStats statsObservable_;
    AutocorrelationStats statsBetterObservable_;
    Statistics statsGeneratorOnBasis_;

    Statistics statsB1_;
    AutocorrelationStats statsB2_;
    Statistics statsD_;
    Vector<double> lastA_;

    AutocorrelationStats statsPostObservable_;
    AutocorrelationStats statsPostBetterObservable_;
    Vector<double> historyObservable_;
    DenseMatrix<double> historyGeneratorOnBasis_;

    Potential* potential_;
  public:
    ControlVariate(Input const& input, Potential& potential, int nbOfFunctions);

    // ACCESSEURS
    virtual bool isNone() const;
    double potential(Vector<double> const& position) const;
    Vector<double> potentialDerivative(Vector<double> const& position) const;
    double potentialLaplacian(Vector<double> const& position) const;
    int decorrelationNbOfIterations() const;
    double decorrelationTime() const;
    virtual int nbOfFunctions() const;
    virtual int nbOfFunctionPairs() const;
		bool doOutput(int iOfIteration) const;
		int const& nbOfAutocoPts() const;
		virtual int nbOfFourier() const;
		virtual int nbOfHermite() const;

    virtual Vector<double> lastB1() const;
    virtual Vector<double> meanB1() const;
    virtual Vector<double> meanB() const;
    virtual DenseMatrix<double> lastD() const;
    virtual DenseMatrix<double> meanD() const;

    virtual double lastObservable() const;
    virtual double meanObservable() const;
    virtual double stdDeviationObservable() const;
    virtual double lastBetterObservable() const;
    virtual double meanBetterObservable() const;
    virtual double stdDeviationBetterObservable() const;

    virtual Vector<double> lastGeneratorOnBasis() const;
    virtual Vector<double> meanGeneratorOnBasis() const;

    virtual double autocorrelation(int indexDifference) const;
    virtual double autocorrelationB2(int indexDifference, int iOfFunction = 0) const;

    virtual Vector<double> lastA() const;
    virtual double lastA(int iOfFunction = 0) const;

    virtual Vector<double> correlationB2() const;
    virtual double correlationB2(int iOfFunction) const;

    // APPEND

    void appendToObservable(double observable, int iOfIteration);
    void appendToB1(double observable, Vector<double>& basisFunction);
    void appendToB2(double observable, Vector<double>& generatorOnBasisFunction, int iOfIteration);
    void appendToD(Vector<double>& generatorOnBasisFunction, Vector<double>& basisFunction);
    void appendToBetterObservable(double observable, Vector<double>& generatorOnBasisFunction, int iOfIteration);

    virtual void update(double observable, Vector<double>& generatorOnBasisFunction, vector<Particle> const& configuration, int iOfIteration);

    // FUNCTION CARACTERIZATION

    virtual double basisFunction(vector<Particle> const& configuration, int iOfFunction = 0) const = 0;


    virtual double laplacianQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const = 0;
    virtual Vector<double> gradientQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const = 0;
    virtual double laplacianP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const = 0;
    virtual Vector<double> gradientP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const = 0;

    Vector<double> generatorHamiltonian(vector<Particle> const& configuration);
    Vector<double> generatorOverdamped(vector<Particle> const& configuration, double beta);
    Vector<double> generatorLangevin(vector<Particle> const& configuration, double beta, double gamma);
    Vector<double> generatorBoundarylangevin(vector<Particle> const& configuration, double betaLeft, double betaRight, double gamma);

    
    
    virtual void display(std::ofstream& out, double time) const;
    virtual void postTreat(std::ofstream& out, double timeStep);

		// TEMP

		virtual void displayMap(ofstream& /*out*/) const{};
		virtual void displayGradQMap(ofstream& /*out*/) const{};
		virtual void displayGradPMap(ofstream& /*out*/) const{};
  };

    class NoControlVariate : public ControlVariate
  {
  public:
    NoControlVariate(Input const& input, Potential& potential);
    virtual int nbOfFunctions() const;
    virtual int nbOfFunctionPairs() const;
    bool isNone() const;
    Vector<double> lastGeneratorOnBasis() const;
    double basisFunction(vector<Particle> const& configuration, int iOfFunction = 0) const;
    virtual double laplacianQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
    virtual Vector<double> gradientQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
    virtual double laplacianP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
    virtual Vector<double> gradientP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
    void update(double observable, Vector<double>& generatorOnBasisFunction, vector<Particle> const& configuration, int iOfIteration);
    virtual void postTreat(std::ofstream& out, double timeStep);
  };


  class SinusControlVariate : public ControlVariate
  {
  public:
    SinusControlVariate(Input const& input, Potential& potential);
    double basisFunction(vector<Particle> const& configuration, int iOfFunction = 0) const;
    virtual double laplacianQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
    virtual Vector<double> gradientQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
    virtual double laplacianP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
    virtual Vector<double> gradientP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
  };

  class CosControlVariate : public ControlVariate
  {
  public:
    CosControlVariate(Input const& input, Potential& potential);
    double basisFunction(vector<Particle> const& configuration, int iOfFunction = 0) const;
    virtual double laplacianQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
    virtual Vector<double> gradientQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
    virtual double laplacianP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
    virtual Vector<double> gradientP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
  };

  class SinExpControlVariate : public ControlVariate
  {
  public:
    SinExpControlVariate(Input const& input, Potential& potential);
    double basisFunction(vector<Particle> const& configuration, int iOfFunction = 0) const;
    virtual double laplacianQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
    virtual Vector<double> gradientQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
    virtual double laplacianP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
    virtual Vector<double> gradientP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
  };

  class CosExpControlVariate : public ControlVariate
  {
  public:
    CosExpControlVariate(Input const& input, Potential& potential);
    double basisFunction(vector<Particle> const& configuration, int iOfFunction = 0) const;
    virtual double laplacianQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
    virtual Vector<double> gradientQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
    virtual double laplacianP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
    virtual Vector<double> gradientP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
  };

  class LangevinControlVariate : public ControlVariate
  {
  public:
    LangevinControlVariate(Input const& input, Potential& potential);
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
    SumEnergyControlVariate(Input const& input, Potential& potential);
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
    EnergyControlVariate(Input const& input, Potential& potential);
    double basisFunction(vector<Particle> const& configuration, int iOfFunction = 0) const;
    virtual double laplacianQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
    virtual Vector<double> gradientQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
    virtual double laplacianP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
    virtual Vector<double> gradientP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
  };

  class LocalControlVariate : public ControlVariate
  {
  public:
    LocalControlVariate(Input const& input, Potential& potential);
    double basisFunction(vector<Particle> const& configuration, int iOfFunction = 0) const;
    virtual double laplacianQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
    virtual Vector<double> gradientQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
    virtual double laplacianP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
    virtual Vector<double> gradientP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
  };

  class KineticControlVariate : public ControlVariate
  {
  public:
    KineticControlVariate(Input const& input, Potential& potential);
    double basisFunction(vector<Particle> const& configuration, int iOfFunction = 0) const;
    virtual double laplacianQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
    virtual Vector<double> gradientQ(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
    virtual double laplacianP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
    virtual Vector<double> gradientP(vector<Particle> const& configuration, int iOfParticle = 0, int iOfFunction = 0) const;
  };



  class TwoControlVariate : public ControlVariate
  {
  public:
    TwoControlVariate(Input const& input, Potential& potential);
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
    BasisControlVariate(const simol::Input& input, simol::Potential& potential, simol::Galerkin* galerkin);
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
		ExpFourierHermiteControlVariate(const simol::Input& input, simol::Potential& potential, Galerkin* galerkin);
		int nbOfFourier() const;
		int nbOfHermite() const;
		virtual void displayMap(ofstream& out) const;
		virtual void displayGradQMap(ofstream& out) const;
		virtual void displayGradPMap(ofstream& out) const;
	};

}

#endif
