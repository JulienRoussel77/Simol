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
    
    //----------- for autocorrelations -------------
    AutocorrelationStats kinTempProfile_;
    AutocorrelationStats potTempTopProfile_;
    AutocorrelationStats potTempBotProfile_;
    AutocorrelationStats bendistProfile_;
    AutocorrelationStats flowProfile_;
    
    //---------- for control variates ------------
    ControlVariate* velocityCV_;
    ControlVariate* forceCV_;
    ControlVariate* lengthCV_;
    ControlVariate* midFlowCV_;
    ControlVariate* sumFlowCV_;

    Output(Input const& input);
    void setControlVariates(Input& input, Potential& potential, Galerkin* galerkin);
    
    ofstream & outObservables();
    ofstream & outParticles();
    ofstream & outParticlesXMakeMol();
    ofstream & outFinalFlow();
    ofstream & outFinalVelocity();
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
    ofstream & outSumFlowPT();
    ofstream & outProfile();
    ofstream & outFinalProfile();
    
    //-- input parameters useful for output --
    const double& timeStep() const;
    double& timeStep();
    double period() const;
    const int& periodNbOfIterations() const;
    const int& longPeriodNbOfIterations() const;
    bool doOutput(int iOfIteration) const;
    bool doLongOutput(int iOfIteration) const;
    const int& nbOfParticles() const;
    const int& nbOfIterations() const;
    double finalTime() const;
    
    //-- fields to output --
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
    
    //-- parametrization of outputs --
    bool doComputeCorrelations() const;
    const int& nbOfAutocoPts() const;
    double autocoPtsPeriod() const;
    const int& decorrelationNbOfIterations() const;
    int& decorrelationNbOfIterations();
    double decorrelationTime() const;
    
    //-- actual outputing functions --
    void displayObservables(int iOfIteration);
    void displayParticles(vector<Particle> const& configuration, int iOfIteration);
    void finalDisplayAutocorrelations();
    void displayFinalVelocity(double temperature, double externalForce, int nbOfFourier = 0, int nbOfHermite = 0);
    
    //-- for NBody systems --
    void displayParticlesXMakeMol(vector<Particle> const& configuration, int iOfIteration, double domainSize=0);
    
    //------------ profiles for chains --------------
    void appendKinTempProfile(double value, int iOfIteration, int iOfParticle);
    void appendPotTempTopProfile(double value, int iOfIteration, int iOfParticle);
    void appendPotTempBotProfile(double value, int iOfIteration, int iOfParticle);
    void appendBendistProfile(double value, int iOfIteration, int iOfParticle);
    void appendFlowProfile(double value, int iOfIteration, int iOfParticle);
    void writeProfile(ofstream & out_, int iOfIteration);
    void displayChainMomenta(vector<Particle> const& configuration, int iOfIteration);
    void displayChainPositions(vector<Particle> const& configuration, int iOfIteration);
    void displayProfile(int iOfIteration);
    void finalChainDisplay(vector<Particle> const& configuration, Vector<double> const& externalForce);
    void displayFinalFlow(double temperature, double delta_temperature, double tau = nan(""), double xi = 0);

    //------------- pour DPDE ---------------
    void displayObservablesDPDE(vector<Particle> const& configuration, int iOfIteration);
 
    //------------- for Galerkin ----------------------
    ControlVariate& velocityCV();
    ControlVariate& forceCV();
    ControlVariate& lengthCV();
    ControlVariate& midFlowCV();
    ControlVariate& sumFlowCV();
    void displayGeneratorOnBasis(ofstream& out, vector<Particle> const& configuration, ControlVariate& controlVariate, double time);
    void updateControlVariate(vector<Particle> const& configuration);

  protected:
    string outputFolderName_;
    std::shared_ptr<ofstream> outObservables_;
    std::shared_ptr<ofstream> outParticles_;
    std::shared_ptr<ofstream> outCorrelation_;
    
    //-- xmakemol outputs for NBody --
    std::shared_ptr<ofstream> outParticlesXMakeMol_;
    
    //-- average velocity for Isolated --
    std::shared_ptr<ofstream> outFinalVelocity_;
    
    //-- control variate outputs --
    std::shared_ptr<ofstream> outVelocitiesGenerator_;
    std::shared_ptr<ofstream> outVelocitiesCV_;
    std::shared_ptr<ofstream> outForcesCV_;
    std::shared_ptr<ofstream> outLengthsCV_;
    
    //-- for chains --
    std::shared_ptr<ofstream> outFinalFlow_;
    std::shared_ptr<ofstream> outBeam_;
    std::shared_ptr<ofstream> outChainVelocities_;
    std::shared_ptr<ofstream> outMidFlowCV_;    
    std::shared_ptr<ofstream> outMidFlowPT_;
    std::shared_ptr<ofstream> outSumFlowCV_;    
    std::shared_ptr<ofstream> outSumFlowPT_;
    std::shared_ptr<ofstream> outProfile_;
    std::shared_ptr<ofstream> outFinalProfile_;
    
    //-- input parameters useful for output --
    int periodNbOfIterations_, longPeriodNbOfIterations_;
    double timeStep_;
    int dimension_;
    int nbOfParticles_;
    int nbOfIterations_;
    double latticeParameter_;
    
    //-- fields to output --
    double kineticEnergy_;
    double potentialEnergy_;
    double internalEnergy_;
    double totalVirial_;
    double energyMidFlow_;
    double energySumFlow_;
    
    //-- parametrization of outputs --
    int decorrelationNbOfIterations_;
    int nbOfAutocoPts_;
    bool doFinalFlow_, doFinalVelocity_;
  };

}
#endif
