#ifndef SIMOL_CONTROLVARIATE_HPP
#define SIMOL_CONTROLVARIATE_HPP

#include "tools.hpp"
#include "statistics.hpp"
#include "potential.hpp"
#include "particle.hpp"
#include "basis.hpp"
#include "galerkin.hpp"


namespace simol
{

  class ControlVariate;
  ControlVariate* createControlVariate(const simol::Input& input, simol::Potential* potential, Galerkin* galerkin, std::size_t iOfReplica = 1);

  class ControlVariate
  {
  friend ControlVariate* createControlVariate(Input const& input, Potential* potential, size_t iOfReplica);
  protected:
    size_t decorrelationNbOfIterations_;
    double decorrelationTime_;
    size_t nbOfFunctions_;
    size_t nbOfFunctionPairs_;
		size_t periodNbOfIterations_;
		int nbOfAutocoPts_;

    AutocorrelationStats<double> statsObservable_;
    AutocorrelationStats<double> statsBetterObservable_;
    Statistics<double> statsGeneratorOnBasis_;

    Statistics<double> statsB1_;
    AutocorrelationStats<double> statsB2_;
    Statistics<double> statsD_;
    VectorXd lastA_;

    AutocorrelationStats<double> statsPostObservable_;
    AutocorrelationStats<double> statsPostBetterObservable_;
    VectorXd historyObservable_;
    MatrixXd historyGeneratorOnBasis_;

    Potential* potential_;
  public:
    ControlVariate(Input const& input, Potential* potential, size_t iOfReplica, size_t nbOfFunctions);

    // ACCESSEURS
    virtual bool isNone() const;
    double potential(Vector<double> const& position) const;
    Vector<double> potentialDerivative(Vector<double> const& position) const;
    double potentialLaplacian(Vector<double> const& position) const;

    size_t decorrelationNbOfIterations() const;
    double decorrelationTime() const;
    virtual size_t nbOfFunctions() const;
    virtual size_t nbOfFunctionPairs() const;
		bool doOutput(size_t iOfIteration) const;
		int const& nbOfAutocoPts() const;
		virtual int nbOfFourier() const;
		virtual int nbOfHermite() const;

    virtual VectorXd lastB1() const;
    virtual VectorXd meanB1() const;
    virtual VectorXd meanB() const;
    virtual MatrixXd lastD() const;
    virtual MatrixXd meanD() const;

    virtual double lastObservable() const;
    virtual double meanObservable() const;
    virtual double stdDeviationObservable() const;
    virtual double lastBetterObservable() const;
    virtual double meanBetterObservable() const;
    virtual double stdDeviationBetterObservable() const;

    virtual VectorXd lastGeneratorOnBasis() const;
    virtual VectorXd meanGeneratorOnBasis() const;

    virtual double autocorrelation(size_t indexDifference) const;
    virtual double autocorrelationB2(size_t indexDifference, size_t iOfFunction = 0) const;

    virtual VectorXd lastA() const;
    virtual double lastA(size_t iOfFunction = 0) const;

    virtual VectorXd correlationB2() const;
    virtual double correlationB2(size_t iOfFunction) const;

    // APPEND

    void appendToObservable(double observable, size_t iOfIteration);
    void appendToB1(double observable, VectorXd& basisFunction);
    void appendToB2(double observable, VectorXd& generatorOnBasisFunction, size_t iOfIteration);
    void appendToD(VectorXd& generatorOnBasisFunction, VectorXd& basisFunction);
    void appendToBetterObservable(double observable, VectorXd& generatorOnBasisFunction, size_t iOfIteration);

    virtual void update(double observable, VectorXd& generatorOnBasisFunction, vector<Particle> const& configuration, size_t iOfIteration);

    // FUNCTION CARACTERIZATION

    virtual double basisFunction(vector<Particle> const& configuration, size_t iOfFunction = 0) const = 0;


    virtual double laplacianQ(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const = 0;
    virtual Vector<double> gradientQ(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const = 0;
    virtual double laplacianP(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const = 0;
    virtual Vector<double> gradientP(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const = 0;

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
    NoControlVariate(Input const& input, Potential* potential, size_t iOfReplica);
    virtual size_t nbOfFunctions() const;
    virtual size_t nbOfFunctionPairs() const;
    bool isNone() const;
    double basisFunction(vector<Particle> const& configuration, size_t iOfFunction = 0) const;
    virtual double laplacianQ(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
    virtual Vector<double> gradientQ(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
    virtual double laplacianP(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
    virtual Vector<double> gradientP(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
    void update(double observable, VectorXd& generatorOnBasisFunction, vector<Particle> const& configuration, size_t iOfIteration);
    virtual void postTreat(std::ofstream& out, double timeStep);
  };


  class SinusControlVariate : public ControlVariate
  {
  public:
    SinusControlVariate(Input const& input, Potential* potential, size_t iOfReplica);
    double basisFunction(vector<Particle> const& configuration, size_t iOfFunction = 0) const;
    virtual double laplacianQ(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
    virtual Vector<double> gradientQ(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
    virtual double laplacianP(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
    virtual Vector<double> gradientP(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
  };

  class CosControlVariate : public ControlVariate
  {
  public:
    CosControlVariate(Input const& input, Potential* potential, size_t iOfReplica);
    double basisFunction(vector<Particle> const& configuration, size_t iOfFunction = 0) const;
    virtual double laplacianQ(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
    virtual Vector<double> gradientQ(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
    virtual double laplacianP(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
    virtual Vector<double> gradientP(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
  };

  class SinExpControlVariate : public ControlVariate
  {
  public:
    SinExpControlVariate(Input const& input, Potential* potential, size_t iOfReplica);
    double basisFunction(vector<Particle> const& configuration, size_t iOfFunction = 0) const;
    virtual double laplacianQ(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
    virtual Vector<double> gradientQ(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
    virtual double laplacianP(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
    virtual Vector<double> gradientP(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
  };

  class CosExpControlVariate : public ControlVariate
  {
  public:
    CosExpControlVariate(Input const& input, Potential* potential, size_t iOfReplica);
    double basisFunction(vector<Particle> const& configuration, size_t iOfFunction = 0) const;
    virtual double laplacianQ(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
    virtual Vector<double> gradientQ(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
    virtual double laplacianP(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
    virtual Vector<double> gradientP(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
  };

  class LangevinControlVariate : public ControlVariate
  {
  public:
    LangevinControlVariate(Input const& input, Potential* potential, size_t iOfReplica);
    double basisFunction(vector<Particle> const& configuration, size_t iOfFunction = 0) const;
    virtual double laplacianQ(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
    virtual Vector<double> gradientQ(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
    virtual double laplacianP(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
    virtual Vector<double> gradientP(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
  };

  class SumEnergyControlVariate : public ControlVariate
  {
    double i0_;
  public:
    SumEnergyControlVariate(Input const& input, Potential* potential, size_t iOfReplica);
    double basisFunction(vector<Particle> const& configuration, size_t iOfFunction = 0) const;
    virtual double laplacianQ(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
    virtual Vector<double> gradientQ(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
    virtual double laplacianP(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
    virtual Vector<double> gradientP(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
  };

  class EnergyControlVariate : public ControlVariate
  {
    double i0_;
  public:
    EnergyControlVariate(Input const& input, Potential* potential, size_t iOfReplica);
    double basisFunction(vector<Particle> const& configuration, size_t iOfFunction = 0) const;
    virtual double laplacianQ(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
    virtual Vector<double> gradientQ(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
    virtual double laplacianP(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
    virtual Vector<double> gradientP(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
  };

  class LocalControlVariate : public ControlVariate
  {
  public:
    LocalControlVariate(Input const& input, Potential* potential, size_t iOfReplica);
    double basisFunction(vector<Particle> const& configuration, size_t iOfFunction = 0) const;
    virtual double laplacianQ(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
    virtual Vector<double> gradientQ(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
    virtual double laplacianP(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
    virtual Vector<double> gradientP(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
  };

  class KineticControlVariate : public ControlVariate
  {
  public:
    KineticControlVariate(Input const& input, Potential* potential, size_t iOfReplica);
    double basisFunction(vector<Particle> const& configuration, size_t iOfFunction = 0) const;
    virtual double laplacianQ(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
    virtual Vector<double> gradientQ(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
    virtual double laplacianP(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
    virtual Vector<double> gradientP(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
  };



  class TwoControlVariate : public ControlVariate
  {
  public:
    TwoControlVariate(Input const& input, Potential* potential, size_t iOfReplica);
    double basisFunction(vector<Particle> const& configuration, size_t iOfFunction = 0) const;
    virtual double laplacianQ(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
    virtual Vector<double> gradientQ(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
    virtual double laplacianP(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
    virtual Vector<double> gradientP(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
  };

	class BasisControlVariate : public ControlVariate
	{
	protected:
		SMat coeffsVec_;
		TensorBasis* basis_;
	public:
    BasisControlVariate(const simol::Input& input, simol::Potential* potential, simol::Galerkin* galerkin, std::size_t iOfReplica);
    double basisFunction(vector<Particle> const& configuration, size_t iOfFunction = 0) const;
    virtual double laplacianQ(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
    virtual Vector<double> gradientQ(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
    virtual double laplacianP(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
    virtual Vector<double> gradientP(vector<Particle> const& configuration, size_t iOfParticle = 0, size_t iOfFunction = 0) const;
	};

	class ExpFourierHermiteControlVariate : public BasisControlVariate
	{
		int nbQ_, nbP_;
		double pMax_, deltaQ_, deltaP_;
	public:
		ExpFourierHermiteControlVariate(const simol::Input& input, simol::Potential* potential, Galerkin* galerkin, std::size_t iOfReplica);
		int nbOfFourier() const;
		int nbOfHermite() const;
		virtual void displayMap(ofstream& out) const;
		virtual void displayGradQMap(ofstream& out) const;
		virtual void displayGradPMap(ofstream& out) const;
	};

}

#endif
