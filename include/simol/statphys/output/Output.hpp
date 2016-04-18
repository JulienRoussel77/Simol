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
      double printPeriodTime() const;
      const int& printPeriodNbOfSteps() const;
      double printLongPeriodTime() const;
      const int& printLongPeriodNbOfSteps() const;
      bool doOutput(int iOfStep) const;
      bool doLongPeriodOutput(int iOfStep) const;
      const int& nbOfParticles() const;
      const int& nbOfSteps() const;
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
      const int& decorrelationNbOfSteps() const;
      int& decorrelationNbOfSteps();
      double decorrelationTime() const;

      //-- actual outputing functions --
      void displayObservables(int iOfStep);
      void displayParticles(vector<Particle> const& configuration, int iOfStep);
      void finalDisplayAutocorrelations();
      void displayFinalVelocity(double temperature, double externalForce, int nbOfFourier = 0, int nbOfHermite = 0);

      //-- for NBody systems --
      void displayParticlesXMakeMol(vector<Particle> const& configuration, int iOfStep, double domainSize = 0);

      //------------ profiles for chains --------------
      void appendKinTempProfile(double value, int iOfStep, int iOfParticle);
      void appendPotTempTopProfile(double value, int iOfStep, int iOfParticle);
      void appendPotTempBotProfile(double value, int iOfStep, int iOfParticle);
      void appendBendistProfile(double value, int iOfStep, int iOfParticle);
      void appendFlowProfile(double value, int iOfStep, int iOfParticle);
      void writeProfile(ofstream & out_, int iOfStep);
      void displayChainMomenta(vector<Particle> const& configuration, int iOfStep);
      void displayChainPositions(vector<Particle> const& configuration, int iOfStep);
      void displayProfile(int iOfStep);
      void finalChainDisplay(vector<Particle> const& configuration, Vector<double> const& externalForce);
      void displayFinalFlow(double temperature, double delta_temperature, double tau = nan(""), double xi = 0);

      //------------- pour DPDE ---------------
      void displayObservablesDPDE(vector<Particle> const& configuration, int iOfStep);

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
      int printPeriodNbOfSteps_, printLongPeriodNbOfSteps_;
      double timeStep_;
      int dimension_;
      int nbOfParticles_;
      int nbOfSteps_;
      double latticeParameter_;

      //-- fields to output --
      double kineticEnergy_;
      double potentialEnergy_;
      double internalEnergy_;
      double totalVirial_;
      double energyMidFlow_;
      double energySumFlow_;

      //-- parametrization of outputs --
      int decorrelationNbOfSteps_;
      int nbOfAutocoPts_;
      bool doFinalFlow_, doFinalVelocity_;
    public:

      //---------- for control variates ------------
      ControlVariate* velocityCV_;
      ControlVariate* forceCV_;
      ControlVariate* lengthCV_;
      ControlVariate* midFlowCV_;
      ControlVariate* sumFlowCV_;

      //----------- for autocorrelations -------------
      AutocorrelationStats kinTempProfile_;
      AutocorrelationStats potTempTopProfile_;
      AutocorrelationStats potTempBotProfile_;
      AutocorrelationStats bendistProfile_;
      AutocorrelationStats flowProfile_;
  };

}
#endif
