#ifndef SIMOL_OUTPUT_HPP
#define SIMOL_OUTPUT_HPP

#include <iomanip>
using std::setw;
#include "Tools.hpp"
#include "Particle.hpp"
#include "Statistics.hpp"
#include "ControlVariate.hpp"
#include "Galerkin.hpp"

namespace simol
{
  
  class Output
  {
  public:
    string outputFolderName_;
    ofstream outObservables_;
    ofstream outParticles_;
    ofstream outParticlesXMakeMol_;
    ofstream outFinalFlow_;
    ofstream outFinalVelocity_;
    ofstream outCorrelation_;
    ofstream outVelocities_;
    ofstream outBeam_;
    ofstream outVelocitiesCV_;
    ofstream outVelocitiesGenerator_;
    ofstream outForcesCV_;
    ofstream outLengthsCV_;
    ofstream outMidFlowCV_;    
    ofstream outMidFlowPT_;
    ofstream outSumFlowCV_;    
    ofstream outSumFlowPT_;
    ofstream outProfile_;
    ofstream outFinalProfile_;
    
    string profilePath_;
    
    size_t periodNbOfIterations_, profilePeriodNbOfIterations_;
    double timeStep_;
    int dimension_;
    size_t nbOfParticles_;
    size_t nbOfIterations_;
    double latticeParameter_;
    
    double kineticEnergy_;
    double potentialEnergy_;
    double internalEnergy_;
    double totalVirial_;
    double energyMidFlow_;
    double energySumFlow_;
    
    size_t decorrelationNbOfIterations_;
    int nbOfAutocoPts_;
    bool doFinalFlow_, doFinalVelocity_;
  public:
    
    ControlVariate* velocityCV_;
    ControlVariate* forceCV_;
    ControlVariate* lengthCV_;
    ControlVariate* midFlowCV_;
    ControlVariate* sumFlowCV_;
    
    AutocorrelationStats kinTempProfile_;
		AutocorrelationStats potTempTopProfile_;
    AutocorrelationStats potTempBotProfile_;
    AutocorrelationStats bendistProfile_;
    AutocorrelationStats flowProfile_;
    
    Output(Input const& input);
    void setControlVariates(Input& input, Potential& potential, Galerkin* galerkin);
    
    const double& timeStep() const;
    double& timeStep();
    double period() const;
    const size_t& periodNbOfIterations() const;
    const size_t& profilePeriodNbOfIterations() const;
    bool doOutput(size_t iOfIteration) const;
    bool doProfileOutput(size_t iOfIteration) const;
    const size_t& nbOfParticles() const;
    const size_t& nbOfIterations() const;
    double finalTime() const;
    
    const double& kineticEnergy() const;
    double& kineticEnergy();
    const double& potentialEnergy() const;
    double& potentialEnergy();
    const double& internalEnergy() const;
    double& internalEnergy();
    const double& totalVirial() const;
    double& totalVirial();
    double energy() const;
    double temperature() const;
    double pressure() const;
    const double& energyMidFlow() const;
    double& energyMidFlow();
    const double& energySumFlow() const;
    double& energySumFlow();
    
    bool doComputeCorrelations() const;
    
    const int& nbOfAutocoPts() const;
    double autocoPtsPeriod() const;
    
    
    /*Vector<double>& responseForces();
      Vector<double> const& responseForces() const;
      double& responseForces(const int& i);
      const double& responseForces(const int& i) const;*/
    const size_t& decorrelationNbOfIterations() const;
    size_t& decorrelationNbOfIterations();
    double decorrelationTime() const;
    
    void writeProfile(ofstream& out_, size_t iOfIteration);
    //double& integratedAutocorrelationP();
    //void display(vector<Particle> const& configuration, size_t iOfIteration);
    void displayObservables(size_t iOfIteration);
    void displayChainMomenta(vector<Particle> const& configuration, size_t iOfIteration);
    void displayChainPositions(vector<Particle> const& configuration, size_t iOfIteration);
    void displayParticles(vector<Particle> const& configuration, size_t iOfIteration);
    void displayParticlesXMakeMol(vector<Particle> const& configuration, size_t iOfIteration, double domainSize=0);
    void displayProfile(size_t iOfIteration);
    
    void finalDisplayAutocorrelations();
    void finalDisplay(vector<Particle> const& configuration, Vector<double> const& externalForce);
    void displayFinalFlow(double temperature, double delta_temperature, double tau = nan(""), double xi = 0);
    void displayFinalVelocity(double temperature, double externalForce, int nbOfFourier = 0, int nbOfHermite = 0);

    //------------- pour Galerkin ----------------------
    ControlVariate& velocityCV();
    ControlVariate& forceCV();
    ControlVariate& lengthCV();
    ControlVariate& midFlowCV();
    ControlVariate& sumFlowCV();
    void displayGeneratorOnBasis(ofstream& out, vector<Particle> const& configuration, ControlVariate& controlVariate, double time);
    void updateControlVariate(vector<Particle> const& configuration);
    
    //------------ profils pour chaines --------------
    void appendKinTempProfile(double value, size_t iOfIteration, size_t iOfParticle);
    void appendPotTempTopProfile(double value, size_t iOfIteration, size_t iOfParticle);
    void appendPotTempBotProfile(double value, size_t iOfIteration, size_t iOfParticle);
    void appendBendistProfile(double value, size_t iOfIteration, size_t iOfParticle);
    void appendFlowProfile(double value, size_t iOfIteration, size_t iOfParticle);
  
    //------------- pour DPDE ---------------
    void displayObservablesDPDE(vector<Particle> const& configuration, size_t iOfIteration);
  };

}
#endif
