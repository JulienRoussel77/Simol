#ifndef SIMOL_OUTPUT_HPP
#define SIMOL_OUTPUT_HPP

#include <iomanip>
using std::setw;
#include "simol/statphys/Tools.hpp"
#include "simol/statphys/system/Particle.hpp"
#include "Statistics.hpp"
#include "simol/statphys/controlVariate/ControlVariate.hpp"
#include "simol/statphys/controlVariate/Galerkin.hpp"
#include "simol/statphys/output/Observable.hpp"
#include "simol/statphys/controlVariate/CVBasis.hpp"

namespace simol
{
  


  class Output
  {

  public:

    Output(Input const& input);
    Observable* addObservable(const Input& input, const string& outPath);
    Observable* addControlVariate(const Input& input, const string& outPath, Galerkin* galerkin);
    void setControlVariates(Input& input, Potential& potential, Galerkin* galerkin);
    

    ofstream & outThermo();
    ofstream & outParticles();
    ofstream & outParticlesXMakeMol();
    ofstream & outParticlesFullConfiguration();
    ofstream & outFinalFlow();
    ofstream & outFinalVelocity();
    ofstream & meanValueObservables();
    ofstream & outCorrelation();
    ofstream & outChainVelocities();
    ofstream & outBeam();
    ofstream & outVelocitiesGenerator();
    ofstream & outProfile();
    ofstream & outFinalProfile();

    //-- input parameters useful for output --
    const double& timeStep() const;
    double& timeStep();
    double printPeriodTime() const;
    int const& printPeriodNbOfSteps() const;
    double printLongPeriodTime() const;
    int const& printLongPeriodNbOfSteps() const;
    bool doOutput(long int iOfStep) const;
    bool doLongPeriodOutput(long int iOfStep) const;
    int const& nbOfParticles() const;
    const long int& nbOfSteps() const;
    double finalTime() const;
    int const& dimension() const;
    double const& latticeParameter() const;
    int nbOfObservables() const;
    vector<Observable*>& observables();
    Observable* observables(int iOfObservable);


    //-- fields to output --
    const double& kineticEnergy() const;
    double& kineticEnergy();
    const double& potentialEnergy() const;
    double& potentialEnergy();
    const double& internalEnergy() const;
    double& internalEnergy();
    const double& internalTemperature() const;
    double& internalTemperature();
    const double& Pressure() const;  
    double& Pressure();
    const double& totalVirial() const;
    double& totalVirial();
    double energy() const;
    double temperature() const;
    double pressure() const;
    const double& rejectionCount() const;
    double& rejectionCount();
    const double& negativeEnergiesCount() const;
    double& negativeEnergiesCount();

    //-- parametrization of outputs --
    bool doComputeCorrelations() const;
    int const& nbOfAutocoPts() const;
    double autocoPtsPeriod() const;
    int const& decorrelationNbOfSteps() const;
    int& decorrelationNbOfSteps();
    double decorrelationTime() const;

    //-- actual outputing functions --
    void displayThermoVariables(long int iOfStep);
    void displayParticles(vector<Particle> const& configuration, long int iOfStep);
    void finalDisplayCorrelations();
    void displayFinalVelocity(double temperature, double externalForce, int nbOfFourier = 0, int nbOfHermite = 0);

    //-- for NBody systems --
    void displayParticlesXMakeMol(vector<Particle> const& configuration, long int iOfStep, double domainSize = 0);
    void displayParticlesFullConfiguration(vector<Particle> const& configuration, long int iOfStep); 

    //------------ profiles for chains --------------
    void appendKinTempProfile(double value, long int iOfStep, int iOfParticle);
    void appendPotTempTopProfile(double value, long int iOfStep, int iOfParticle);
    void appendPotTempBotProfile(double value, long int iOfStep, int iOfParticle);
    void appendBendistProfile(double value, long int iOfStep, int iOfParticle);
    void appendFlowProfile(double value, long int iOfStep, int iOfParticle);
    void writeProfile(ofstream & out_, long int iOfStep);
    void displayChainMomenta(vector<Particle> const& configuration, long int iOfStep);
    void displayChainPositions(vector<Particle> const& configuration, long int iOfStep);
    void displayProfile(long int iOfStep);
    void finalChainDisplay(vector<Particle> const& configuration, Vector<double> const& externalForce);
    void displayFinalFlow(double temperature, double delta_temperature, double parameter1 = 0, double parameter2 = 0);

      //------------- pour DPDE ---------------
    void displayThermoVariablesDPDE(vector<Particle> const& configuration, long int iOfStep);
    void finalDisplayCorrelationsDPDE();
    
    Observable& obsKineticEnergy();
    Observable const& obsKineticEnergy() const;
    Observable& obsPotentialEnergy();
    Observable const& obsPotentialEnergy() const;
    Observable& obsPressure();
    Observable const& obsPressure() const;
    Observable& obsInternalEnergy();
    Observable const& obsInternalEnergy() const;
    Observable& obsInternalTemperature();
    Observable const& obsInternalTemperature() const;
    Observable& obsVelocity();
    Observable const& obsVelocity() const;
    Observable& obsForce();
    Observable const& obsForce() const;
    Observable& obsLength();
    Observable const& obsLength() const;
    Observable& obsMidFlow();
    Observable const& obsMidFlow() const;
    Observable& obsSumFlow();
    Observable const& obsSumFlow() const;
    Observable& obsModiFlow();
    Observable const& obsModiFlow() const;

    bool hasControlVariate() const;
  protected:
    string outputFolderName_;
    std::shared_ptr<ofstream> outThermo_;
    std::shared_ptr<ofstream> outParticles_;
    std::shared_ptr<ofstream> outCorrelation_;

    //-- configuration outputs for NBody --
    std::shared_ptr<ofstream> outParticlesXMakeMol_;
    std::shared_ptr<ofstream> outParticlesFullConfiguration_;

    //-- average velocity for Isolated --
    std::shared_ptr<ofstream> outFinalVelocity_;

    //-- mean observables for DPDE --
    std::shared_ptr<ofstream> meanValueObservables_;

    //-- control variate outputs --
    std::shared_ptr<ofstream> outVelocitiesGenerator_;

    //-- for chains --
    std::shared_ptr<ofstream> outFinalFlow_;
    std::shared_ptr<ofstream> outBeam_;  
    std::shared_ptr<ofstream> outChainVelocities_;
    std::shared_ptr<ofstream> outProfile_;
    std::shared_ptr<ofstream> outFinalProfile_;

    //-- input parameters useful for output --
    int printPeriodNbOfSteps_, printLongPeriodNbOfSteps_;
    double timeStep_;
    int dimension_;
    int nbOfParticles_;
    long int nbOfSteps_;
    double latticeParameter_;

    //-- fields to output --
    double totalVirial_;
    double rejectionCount_;
    double negativeEnergiesCount_;

    //-- parametrization of outputs --
    int decorrelationNbOfSteps_;
    int nbOfAutocoPts_;
    bool doFinalFlow_, doFinalVelocity_;
  public:
    std::shared_ptr<ofstream> outTest_;   // for debug purpose only
    
    Observable* obsKineticEnergy_;
    Observable* obsPotentialEnergy_;
    Observable* obsPressure_;
    Observable* obsInternalEnergy_;
    Observable* obsInternalTemperature_;
    Observable* obsVelocity_;
    Observable* obsForce_;
    Observable* obsLength_;
    Observable* obsMidFlow_;
    Observable* obsSumFlow_;
    Observable* obsModiFlow_;
    
    vector<Observable*> observables_;

    //----------- for autocorrelations -------------
    //-- chains --
    AutocorrelationStats kinTempProfile_;
    AutocorrelationStats potTempTopProfile_;
    AutocorrelationStats potTempBotProfile_;
    AutocorrelationStats bendistProfile_;
    AutocorrelationStats flowProfile_;
    
    CVBasis cvBasis_;
    
  };
  

  
}
#endif
