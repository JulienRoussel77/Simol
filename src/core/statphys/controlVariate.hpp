#pragma once

#include "statistics.hpp"
#include "potential.hpp"
#include "particle.hpp"

 
namespace simol
{
  
  class ControlVariate;
  ControlVariate* createControlVariate(Input const& input, Potential& potential, size_t indexOfReplica=1);
  
  class ControlVariate
  {
  friend ControlVariate* createControlVariate(Input const& input, Potential& potential, size_t indexOfReplica);
  protected:
    size_t decorrelationNumberOfIterations_;
    double decorrelationTime_;
    size_t nbOfFunctions_;
    size_t nbOfFunctionPairs_;
    
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
    ControlVariate(Input const& input, Potential& potential, size_t indexOfReplica, size_t nbOfFunctions);
    
    // ACCESSEURS
    virtual bool isNone() const;
    double potential(dvec const& position) const;
    dvec potentialDerivative(dvec const& position) const;
    double potentialLaplacian(dvec const& position) const;
    
    size_t decorrelationNumberOfIterations() const;
    double decorrelationTime() const;
    virtual size_t nbOfFunctions() const;
    virtual size_t nbOfFunctionPairs() const;
    
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
    
    void appendToObservable(double observable, size_t indexOfIteration);
    void appendToB1(double observable, VectorXd& basisFunction);
    void appendToB2(double observable, VectorXd& generatorOnBasisFunction, size_t indexOfIteration);
    void appendToD(VectorXd& generatorOnBasisFunction, VectorXd& basisFunction);
    void appendToBetterObservable(double observable, VectorXd& generatorOnBasisFunction, size_t indexOfIteration);
   
    virtual void update(double observable, VectorXd& generatorOnBasisFunction, vector<Particle> const& configuration, size_t indexOfIteration);    
    
    // FUNCTION CARACTERIZATION
    
    virtual double basisFunction(vector<Particle> const& configuration, size_t iOfFunction = 0) const = 0;

    
    virtual double laplacienQ(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const = 0;
    virtual dvec gradientQ(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const = 0;
    virtual double laplacienP(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const = 0;
    virtual dvec gradientP(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const = 0;
  
    virtual void display(std::ofstream& out, double time) const;
    virtual void postTreat(std::ofstream& out, double timeStep);
  };
  
    class NoControlVariate : public ControlVariate
  {
  public:
    NoControlVariate(Input const& input, Potential& potential, size_t indexOfReplica);
    virtual size_t nbOfFunctions() const;
    virtual size_t nbOfFunctionPairs() const;
    bool isNone() const;
    double basisFunction(vector<Particle> const& configuration, size_t iOfFunction = 0) const;
    virtual double laplacienQ(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
    virtual dvec gradientQ(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
    virtual double laplacienP(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
    virtual dvec gradientP(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
    void update(double observable, VectorXd& generatorOnBasisFunction, vector<Particle> const& configuration, size_t indexOfIteration);
    virtual void postTreat(std::ofstream& out, double timeStep);
  };
  
  
  class SinusControlVariate : public ControlVariate
  {
  public:
    SinusControlVariate(Input const& input, Potential& potential, size_t indexOfReplica);
    double basisFunction(vector<Particle> const& configuration, size_t iOfFunction = 0) const;
    virtual double laplacienQ(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
    virtual dvec gradientQ(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
    virtual double laplacienP(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
    virtual dvec gradientP(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
  };
    
  class CosControlVariate : public ControlVariate
  {
  public:
    CosControlVariate(Input const& input, Potential& potential, size_t indexOfReplica);
    double basisFunction(vector<Particle> const& configuration, size_t iOfFunction = 0) const;
    virtual double laplacienQ(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
    virtual dvec gradientQ(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
    virtual double laplacienP(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
    virtual dvec gradientP(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
  };
  
  class SinExpControlVariate : public ControlVariate
  {
  public:
    SinExpControlVariate(Input const& input, Potential& potential, size_t indexOfReplica);
    double basisFunction(vector<Particle> const& configuration, size_t iOfFunction = 0) const;
    virtual double laplacienQ(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
    virtual dvec gradientQ(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
    virtual double laplacienP(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
    virtual dvec gradientP(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
  };
  
  class CosExpControlVariate : public ControlVariate
  {
  public:
    CosExpControlVariate(Input const& input, Potential& potential, size_t indexOfReplica);
    double basisFunction(vector<Particle> const& configuration, size_t iOfFunction = 0) const;
    virtual double laplacienQ(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
    virtual dvec gradientQ(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
    virtual double laplacienP(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
    virtual dvec gradientP(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
  };
  
  class LangevinControlVariate : public ControlVariate
  {
  public:
    LangevinControlVariate(Input const& input, Potential& potential, size_t indexOfReplica);
    double basisFunction(vector<Particle> const& configuration, size_t iOfFunction = 0) const;
    virtual double laplacienQ(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
    virtual dvec gradientQ(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
    virtual double laplacienP(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
    virtual dvec gradientP(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
  };
  
  class SumEnergyControlVariate : public ControlVariate
  {
    double i0_;
  public:
    SumEnergyControlVariate(Input const& input, Potential& potential, size_t indexOfReplica);
    double basisFunction(vector<Particle> const& configuration, size_t iOfFunction = 0) const;
    virtual double laplacienQ(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
    virtual dvec gradientQ(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
    virtual double laplacienP(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
    virtual dvec gradientP(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
  };  
  
  class EnergyControlVariate : public ControlVariate
  {
    double i0_;
  public:
    EnergyControlVariate(Input const& input, Potential& potential, size_t indexOfReplica);
    double basisFunction(vector<Particle> const& configuration, size_t iOfFunction = 0) const;
    virtual double laplacienQ(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
    virtual dvec gradientQ(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
    virtual double laplacienP(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
    virtual dvec gradientP(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
  };
  
  class LocalControlVariate : public ControlVariate
  {
  public:
    LocalControlVariate(Input const& input, Potential& potential, size_t indexOfReplica);
    double basisFunction(vector<Particle> const& configuration, size_t iOfFunction = 0) const;
    virtual double laplacienQ(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
    virtual dvec gradientQ(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
    virtual double laplacienP(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
    virtual dvec gradientP(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
  };
  
  class KineticControlVariate : public ControlVariate
  {
  public:
    KineticControlVariate(Input const& input, Potential& potential, size_t indexOfReplica);
    double basisFunction(vector<Particle> const& configuration, size_t iOfFunction = 0) const;
    virtual double laplacienQ(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
    virtual dvec gradientQ(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
    virtual double laplacienP(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
    virtual dvec gradientP(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
  };
  
  
  
  class TwoControlVariate : public ControlVariate
  {
  public:
    TwoControlVariate(Input const& input, Potential& potential, size_t indexOfReplica);
    double basisFunction(vector<Particle> const& configuration, size_t iOfFunction = 0) const;
    virtual double laplacienQ(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
    virtual dvec gradientQ(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
    virtual double laplacienP(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
    virtual dvec gradientP(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
  };
  
  class IOControlVariate : public ControlVariate
  {
  public:
    IOControlVariate(Input const& input, Potential& potential, size_t indexOfReplica);
    double basisFunction(vector<Particle> const& configuration, size_t iOfFunction = 0) const;
    virtual double laplacienQ(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
    virtual dvec gradientQ(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
    virtual double laplacienP(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
    virtual dvec gradientP(vector<Particle> const& configuration, size_t indexOfParticle = 0, size_t iOfFunction = 0) const;
  };
  
}