#ifndef SIMOL_OUTPUT_HPP
#define SIMOL_OUTPUT_HPP


#include "simol/statphys/Tools.hpp"
#include "simol/statphys/system/Particle.hpp"
#include "simol/statphys/system/System.hpp"
#include "Statistics.hpp"
#include "simol/statphys/controlVariate/ControlVariate.hpp"
#include "simol/statphys/controlVariate/Galerkin.hpp"
#include "simol/statphys/output/Observable.hpp"
#include "simol/statphys/controlVariate/CVBasis.hpp"
#include "simol/statphys/dynamics/DynamicsParameters.hpp"

namespace simol
{

  class Output
  {

  public:

    Output(Input const& input, shared_ptr<CVBasis> cvBasis0);
    virtual ~Output();
      
    Observable*& getObservablePtr(int idObs);
    void addObservable(const Input& input, int idObs, int decorrelationNbOfSteps0, int nbOfAutocoPts0);
    

    ofstream & outThermo() {return *outThermo_;}
    ofstream & outParticles() {return *outParticles_;}
    ofstream & outXMakeMol() {return *outXMakeMol_;}
    ofstream & outBackUp() {return *outBackUp_;}
    ofstream & outFinalFlux() {return *outFinalFlux_;}
    ofstream & outFinalLength() {return *outFinalLength_;}
    ofstream & outFinalVelocity() {return *outFinalVelocity_;}
    ofstream & outFinalLagrangeMultiplier() {return *outFinalLagrangeMultiplier_;}
    ofstream & outMeanThermo() {return *outMeanThermo_;}
    //ofstream & outChainVelocities() {return *outChainVelocities_;}
    //ofstream & outBeam() {return *outBeam_;}
    ofstream & outVelocitiesGenerator() {return *outVelocitiesGenerator_;}
    ofstream & outInstantProfile() {return *outInstantProfile_;}
    ofstream & outProfile() {return *outProfile_;}
    ofstream & outFinalProfile() {return *outFinalProfile_;}

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
    
    //-- test if a feature is activated --
    bool doOutParticles() const {return (bool)outParticles_;}    
    bool doDPDE() const {return doDPDE_;} 
    bool doXMakeMol() const {return doXMakeMol_;}    
    bool doOutBackUp() const {return (bool)outBackUp_;}
    bool doOutChain() const {return (bool)doOutChain_;}
    
    bool doFinalLength() const {return doFinalLength_;}
    bool doFinalVelocity() const {return doFinalVelocity_;}
    bool doFinalLagrangeMultiplier() const {return doFinalLagrangeMultiplier_;}
    bool doFinalFlux() const {return doFinalFlux_;}
    
    //-- fields to output --
    const double& kineticEnergy() const;
    double& kineticEnergy();
    const double& potentialEnergy() const;
    double& potentialEnergy();
    const double& internalEnergy() const;
    double& internalEnergy();
    const double& internalTemperature() const;
    double& internalTemperature();
    const double& length() const;  
    double& length();
    const double& velocity() const;  
    double& velocity();
    const double& force() const;  
    double& force();
    const double& pressure() const;  
    double& pressure();
    const double& totalVirial() const;
    double& totalVirial();
    const double& totalEnergy() const;
    double& totalEnergy();
    const double& temperature() const;
    double& temperature();
    const double& lagrangeMultiplier() const;
    double& lagrangeMultiplier();
    //double temperature() const;
    //double pressure() const;
    const double& rejectionCountFD() const;
    double& rejectionCountFD();
    const double& negativeEnergiesCountFD() const;
    double& negativeEnergiesCountFD();
    const double& rejectionCountThermal() const;
    double& rejectionCountThermal();
    const double& negativeEnergiesCountThermal() const;
    double& negativeEnergiesCountThermal();

    //-- parametrization of outputs --
    int const& nbOfAutocoPts() const;
    double autocoPtsPeriod() const;
    int const& nbOfShortAutocoPts() const;
    double shortAutocoPtsPeriod() const;
    int const& decorrelationNbOfSteps() const;
    //int& decorrelationNbOfSteps();
    double decorrelationTime() const;
    int const& shortDecorrelationNbOfSteps() const;
    double shortDecorrelationTime() const;

    //-- actual outputing functions --
    void displayThermoVariables(long int iOfStep);
    void displayParticles(System const& syst, long int iOfStep);
    void finalDisplayCorrelations();
    void displayFinalLength();
    void displayFinalVelocity();

    //-- for NBody systems --
    void displayXMakeMol(System const& syst, long int iOfStep, double domainSize = 0);
    void displayBackUp(System const& syst, long int iOfStep); 

    //------------ profiles for chains --------------
    void appendKinTempProfile(double value, long int iOfStep, int iOfParticle);
    void appendPotTempTopProfile(double value, long int iOfStep, int iOfParticle);
    void appendPotTempBotProfile(double value, long int iOfStep, int iOfParticle);
    void appendBendistProfile(double value, long int iOfStep, int iOfParticle);
    void appendFluxProfile(double value, long int iOfStep, int iOfParticle);
    void appendModiFluxProfile(double value, long int iOfStep, int iOfParticle);
    void appendExtFluxProfile(double value, long int iOfStep, int iOfParticle);
    
    void writeInstantProfile(System const& syst, ofstream & out_, long int iOfStep);
    void writeProfile(ofstream & out_, long int iOfStep);
    //void displayChainMomenta(System const& syst, long int iOfStep);
    //void displayChainPositions(System const& syst, long int iOfStep);
    void displayInstantProfile(System const& syst, long int iOfStep);
    void displayProfile(long int iOfStep);
    void finalChainDisplay();
    void displayFinalFlux(double parameter1=0, double parameter2=0, double parameter3=0);
    void displayFinalChainLagrangeMultiplier(double parameter1=0, double parameter2=0, double parameter3=0);
    
    //-------------- ConstrainedLangevin ----------------
    void displayFinalLagrangeMultiplier();

    //------------- pour DPDE ---------------
    //void finalDisplayCorrelationsDPDE();
    
    Observable& obsKineticEnergy();
    Observable const& obsKineticEnergy() const;
    Observable& obsPotentialEnergy();
    Observable const& obsPotentialEnergy() const;
    Observable& obsTotalEnergy();
    Observable const& obsTotalEnergy() const;
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
    Observable& obsMidFlux();
    Observable const& obsMidFlux() const;
    Observable& obsSumFlux();
    Observable const& obsSumFlux() const;
    Observable& obsModiFlux();
    Observable const& obsModiFlux() const;
    Observable& obsLagrangeMultiplier();
    Observable const& obsLagrangeMultiplier() const;

    bool hasControlVariate() const;
  protected:
    string outputFolderName_;
    
    // output fluxes
    std::shared_ptr<ofstream> outThermo_;
    std::shared_ptr<ofstream> outMeanThermo_;
    std::shared_ptr<ofstream> outParticles_;
    std::shared_ptr<ofstream> outCorrelation_;

    //-- configuration outputs for NBody --
    std::shared_ptr<ofstream> outXMakeMol_;
    std::shared_ptr<ofstream> outBackUp_;

    //-- average velocity for Isolated --
    std::shared_ptr<ofstream> outFinalLength_;
    std::shared_ptr<ofstream> outFinalVelocity_;
    std::shared_ptr<ofstream> outFinalLagrangeMultiplier_;

    //-- control variate outputs --
    std::shared_ptr<ofstream> outVelocitiesGenerator_;

    //-- for chains --
    std::shared_ptr<ofstream> outFinalFlux_;
    std::shared_ptr<ofstream> outInstantProfile_;
    std::shared_ptr<ofstream> outProfile_;
    std::shared_ptr<ofstream> outFinalProfile_;

    // number of steps between two outputs
    // heavy outputs are printed less often than lighter ones
    int printPeriodNbOfSteps_, printLongPeriodNbOfSteps_;
    // time step of the integrator of the dynamics
    double timeStep_;
    // dimension of the system (usually 1, 2 or 3)
    int dimension_;
    // number of particles
    int nbOfParticles_;
    // number of integration iterations
    long int nbOfSteps_;
    // for systems initialized on a lattice, this is the size of one cell
    double latticeParameter_;
    
  public:
    // contains the parameters of the dynamics
    DynamicsParameters parameters_;
    
  protected:

    //-- fields to output --
    double totalVirial_;
    double temperature_;
    // rejection rates (DPDE)
    double rejectionCountFD_;
    double negativeEnergiesCountFD_;
    double rejectionCountThermal_;
    double negativeEnergiesCountThermal_;

    // number of time steps corresponding to the autocorrelatio profile of an observable
    // in practice it should be an upper bound of the decorrelation time in order to get a good estimation of the asymptotic variance
    // some observables have faster decorrelation and they require a shorter profile
    int decorrelationNbOfSteps_, shortDecorrelationNbOfSteps_;
    // number of nodes in the final plot of the correlation profile (number of bins)
    int nbOfAutocoPts_, nbOfShortAutocoPts_;
    // Is true if the corresponding output should be printed
    bool doOutChain_, doFinalFlux_, doFinalLength_, doFinalVelocity_, doFinalLagrangeMultiplier_, doDPDE_, doXMakeMol_;
    // For a thermal chain: is true if the harmonic potential fitting the interaction potential should be fitted numerically (used for control variates, cf paper SIAM:MMS)
    bool fitModifFlux_;
    
  public:
    // for debug purpose only
    std::shared_ptr<ofstream> outTest_;
    // redundent pointers to the observables which are used
    // if all goes well "observables_[idOfObs] == getObservable(idObs)"
    // the observables which are not relevent are not created and remain nullptr
    //vector<Observable*> allObservables_;
    Observable* obsKineticEnergy_;
    Observable* obsPotentialEnergy_;
    Observable* obsTotalEnergy_;
    Observable* obsPressure_;
    Observable* obsInternalEnergy_;
    Observable* obsInternalTemperature_;
    Observable* obsVelocity_;
    Observable* obsForce_;
    Observable* obsLength_;
    Observable* obsMidFlux_;
    Observable* obsSumFlux_;
    Observable* obsModiFlux_;
    Observable* obsLagrangeMultiplier_;
    
    vector<Observable*> observables_;

    //----------- for autocorrelations -------------
    // objects making statistics on the values gathered for the observable profile in a chain
    // allows for example to plot the profile of mean fluxes with an error bar for each site
    AutocorrelationStats kinTempProfile_;
    AutocorrelationStats potTempTopProfile_;
    AutocorrelationStats potTempBotProfile_;
    AutocorrelationStats bendistProfile_;
    AutocorrelationStats fluxProfile_;
    AutocorrelationStats modiFluxProfile_;
    AutocorrelationStats extFluxProfile_;
    
    shared_ptr<CVBasis> cvBasis_;
  };
  

  
}
#endif
