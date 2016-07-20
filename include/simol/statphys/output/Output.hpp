#ifndef SIMOL_OUTPUT_HPP
#define SIMOL_OUTPUT_HPP

#include <iomanip>
using std::setw;
#include "simol/statphys/Tools.hpp"
#include "simol/statphys/system/Particle.hpp"
#include "Statistics.hpp"
#include "simol/statphys/controlVariate/ControlVariate.hpp"
#include "simol/statphys/controlVariate/Galerkin.hpp"

namespace simol
{

  class Output
  {

  public:

    Output(Input const& input);
    void setControlVariates(Input& input, Potential& potential, Galerkin* galerkin);

    ofstream & outObservables();
    ofstream & outParticles();
    ofstream & outParticlesXMakeMol();
    ofstream & outParticlesFullConfiguration();
    ofstream & outFinalFlow();
    ofstream & outFinalVelocity();
    ofstream & meanValueObservables();
    ofstream & outCorrelation();
    ofstream & outChainVelocities();
    ofstream & outBeam();
    ofstream & outVelocitiesCV();
    ofstream & outVelocitiesGenerator();
    ofstream & outForcesCV();
    ofstream & outLengthsCV();
    ofstream & outMidFlowCV();
    ofstream & outMidFlowPT();
    ofstream & outSumFlowCV();
    ofstream & outModiFlowCV();
    ofstream & outSumFlowPT();
    ofstream & outProfile();
    ofstream & outFinalProfile();

    //-- input parameters useful for output --
    const double& timeStep() const;
    double& timeStep();
    double printPeriodTime() const;
    const int& printPeriodNbOfSteps() const;
    double printLongPeriodTime() const;
    const int& printLongPeriodNbOfSteps() const;
    bool doOutput(long int iOfStep) const;
    bool doLongPeriodOutput(long int iOfStep) const;
    const int& nbOfParticles() const;
    const long int& nbOfSteps() const;
    double finalTime() const;

    //-- fields to output --
    const double& kineticEnergy() const;
    double& kineticEnergy();
    const double& potentialEnergy() const;
    double& potentialEnergy();
    const double& internalEnergy() const;
    double& internalEnergy();
    const double& internalTemperature() const;
    double& internalTemperature();
    const double& totalVirial() const;
    double& totalVirial();
    double energy() const;
    double temperature() const;
    double pressure() const;
    const double& energyMidFlow() const;
    double& energyMidFlow();
    const double& energySumFlow() const;
    double& energySumFlow();
    const double& energyModiFlow() const;
    double& energyModiFlow();
    const double& rejectionCount() const;
    double& rejectionCount();

    //-- parametrization of outputs --
    bool doComputeCorrelations() const;
    const int& nbOfAutocoPts() const;
    double autocoPtsPeriod() const;
    const int& decorrelationNbOfSteps() const;
    int& decorrelationNbOfSteps();
    double decorrelationTime() const;

    //-- actual outputing functions --
    void displayObservables(long int iOfStep);
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
    void displayObservablesDPDE(vector<Particle> const& configuration, long int iOfStep);
    void appendKineticEnergy(double value, long int iOfStep);
    void appendPotentialEnergy(double value, long int iOfStep);
    void appendInternalEnergy(double value, long int iOfStep);
    void appendInternalTemperature(double value, long int iOfStep);
    void appendPressure(double value, long int iOfStep);
    void finalDisplayCorrelationsDPDE();

    //------------- for Galerkin ----------------------
    ControlVariate& velocityCV();
    ControlVariate& forceCV();
    ControlVariate& lengthCV();
    ControlVariate& midFlowCV();
    ControlVariate& sumFlowCV();
    ControlVariate& modiFlowCV();
    void displayGeneratorOnBasis(ofstream& out, vector<Particle> const& configuration, ControlVariate& controlVariate, double time);
    void updateControlVariate(vector<Particle> const& configuration);

  protected:
    string outputFolderName_;
    std::shared_ptr<ofstream> outObservables_;
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
    std::shared_ptr<ofstream> outVelocitiesCV_;
    std::shared_ptr<ofstream> outForcesCV_;

    //-- for chains --
    std::shared_ptr<ofstream> outFinalFlow_;
    std::shared_ptr<ofstream> outBeam_;
    std::shared_ptr<ofstream> outLengthsCV_;      
    std::shared_ptr<ofstream> outChainVelocities_;
    std::shared_ptr<ofstream> outMidFlowCV_;
    std::shared_ptr<ofstream> outMidFlowPT_;
    std::shared_ptr<ofstream> outSumFlowCV_;
    std::shared_ptr<ofstream> outSumFlowPT_;
    std::shared_ptr<ofstream> outModiFlowCV_;
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
    double kineticEnergy_;
    double potentialEnergy_;
    double internalEnergy_;
    double totalVirial_;
    double energyMidFlow_;
    double energySumFlow_;
    double energyModiFlow_;
    double rejectionCount_;
    double internalTemperature_;

    //-- parametrization of outputs --
    int decorrelationNbOfSteps_;
    int nbOfAutocoPts_;
    bool doFinalFlow_, doFinalVelocity_;
  public:
    std::shared_ptr<ofstream> outTest_;   // for debug purpose only

    //---------- for control variates ------------
    ControlVariate* velocityCV_;
    ControlVariate* forceCV_;
    ControlVariate* lengthCV_;
    ControlVariate* midFlowCV_;
    ControlVariate* sumFlowCV_;
    ControlVariate* modiFlowCV_;

    //----------- for autocorrelations -------------
    //-- chains --
    AutocorrelationStats kinTempProfile_;
    AutocorrelationStats potTempTopProfile_;
    AutocorrelationStats potTempBotProfile_;
    AutocorrelationStats bendistProfile_;
    AutocorrelationStats flowProfile_;
    
    //-- DPDE --
    AutocorrelationStats averageKineticEnergy_;
    AutocorrelationStats averagePotentialEnergy_;
    AutocorrelationStats averageInternalEnergy_;
    AutocorrelationStats averageInternalTemperature_;
    AutocorrelationStats averagePressure_;
    
  };
  
}
#endif
