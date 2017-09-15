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

namespace simol
{

  class Output
  {

  public:

    Output(Input const& input, shared_ptr<CVBasis> cvBasis0);
    virtual ~Output();
      
    //Observable* addObservable(const Input& input, const string& outPath);
    Observable*& getObservablePtr(int idObs);
    void addObservable(const Input& input, int idObs, int decorrelationNbOfSteps0, int nbOfAutocoPts0);
    //Observable* addControlVariate(const Input& input, const string& outPath, Galerkin* galerkin);
    //void setControlVariates(Input& input, Potential& potential, Galerkin* galerkin);
    

    ofstream & outThermo() {return *outThermo_;}
    ofstream & outParticles() {return *outParticles_;}
    ofstream & outXMakeMol() {return *outXMakeMol_;}
    ofstream & outBackUp() {return *outBackUp_;}
    ofstream & outFinalFlow() {return *outFinalFlow_;}
    ofstream & outFinalLength() {return *outFinalLength_;}
    ofstream & outFinalVelocity() {return *outFinalVelocity_;}
    ofstream & outFinalLagrangeMultiplier() {return *outFinalLagrangeMultiplier_;}
    ofstream & outMeanThermo() {return *outMeanThermo_;}
    ofstream & outChainVelocities() {return *outChainVelocities_;}
    ofstream & outBeam() {return *outBeam_;}
    ofstream & outVelocitiesGenerator() {return *outVelocitiesGenerator_;}
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
    bool doOutChain() const {return (bool)outProfile_;}
    
    /*bool doFinalLength() const {return (bool)outFinalLength_;}
    bool doFinalVelocity() const {return (bool)outFinalVelocity_;}
    bool doFinalLagrangeMultiplier() const {return (bool)outFinalLagrangeMultiplier_;}
    bool doFinalChainLagrangeMultiplier() const {return (bool)outFinalChainLagrangeMultiplier_ && ;}
    bool doFinalFlow() const {return (bool)outFinalFlow_;}*/
    
    bool doFinalLength() const {return doFinalLength_;}
    bool doFinalVelocity() const {return doFinalVelocity_;}
    bool doFinalLagrangeMultiplier() const {return doFinalLagrangeMultiplier_;}
    bool doFinalChainLagrangeMultiplier() const {return doFinalChainLagrangeMultiplier_;}
    bool doFinalFlow() const {return doFinalFlow_;}
    
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
    void appendFlowProfile(double value, long int iOfStep, int iOfParticle);
    void appendModiFlowProfile(double value, long int iOfStep, int iOfParticle);
    void writeProfile(ofstream & out_, long int iOfStep);
    void displayChainMomenta(System const& syst, long int iOfStep);
    void displayChainPositions(System const& syst, long int iOfStep);
    void displayProfile(long int iOfStep);
    void finalChainDisplay();
    void displayFinalFlow(double parameter1=0, double parameter2=0, double parameter3=0);
    void displayFinalChainLagrangeMultiplier(double parameter1=0, double parameter2=0, double parameter3=0);
    
    //-------------- ConstrainedLangevin ----------------
    void displayFinalLagrangeMultiplier();

    //------------- pour DPDE ---------------
    //void finalDisplayCorrelationsDPDE();
    
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
    Observable& obsLagrangeMultiplier();
    Observable const& obsLagrangeMultiplier() const;

    bool hasControlVariate() const;
  protected:
    string outputFolderName_;
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
    
  public:
    double constTemperature_;
    double constTemperatureLeft_;
    double constTemperatureRight_;
    double constDeltaTemperature_;
    double constGamma_;
    double constXi_;
    double constTauBending_;
    double constNonEqAmplitude_;
    double constInteractionRatio_;
    int constNbOfQModes_;
    int constNbOfPModes_;
    double constDrift_;
    double const constBulkDriving_;
    double const constFlux_;
    
  protected:

    //-- fields to output --
    double totalEnergy_;
    double totalVirial_;
    double temperature_;
    // rejection rates (DPDE)
    double rejectionCountFD_;
    double negativeEnergiesCountFD_;
    double rejectionCountThermal_;
    double negativeEnergiesCountThermal_;

    //-- parametrization of outputs --
    int decorrelationNbOfSteps_, shortDecorrelationNbOfSteps_;
    
    int nbOfAutocoPts_, nbOfShortAutocoPts_;
    bool doFinalFlow_, doFinalLength_, doFinalVelocity_, doFinalLagrangeMultiplier_, doFinalChainLagrangeMultiplier_, doDPDE_, doXMakeMol_;
    bool fitModifFlow_;
    
  public:
    
    std::shared_ptr<ofstream> outTest_;   // for debug purpose only
    
    vector<Observable*> allObservables_;
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
    Observable* obsLagrangeMultiplier_;
    
    vector<Observable*> observables_;

    //----------- for autocorrelations -------------
    //-- chains --
    AutocorrelationStats kinTempProfile_;
    AutocorrelationStats potTempTopProfile_;
    AutocorrelationStats potTempBotProfile_;
    AutocorrelationStats bendistProfile_;
    AutocorrelationStats flowProfile_;
    AutocorrelationStats modiFlowProfile_;
    
    shared_ptr<CVBasis> cvBasis_;
  };
  

  
}
#endif
