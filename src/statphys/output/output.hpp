#ifndef SIMOL_OUTPUT_HPP
#define SIMOL_OUTPUT_HPP

#include <iomanip>
using std::setw;
#include "tools.hpp"
#include "particle.hpp"
#include "statistics.hpp"
#include "controlVariate.hpp"
#include "galerkin.hpp"

namespace simol
{
  


  
  class Output
  {
	public:
    string outputFolderName_;
    ofstream outObservables_;
    ofstream outParticles_;
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
    
    int verbose_;
    size_t periodNbOfIterations_, profilePeriodNbOfIterations_;
    double timeStep_;
    int dimension_;
    size_t nbOfParticles_;
		size_t nbOfIterations_;
    
    double kineticEnergy_;
    double potentialEnergy_;
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
		
		AutocorrelationStats<double> kinTempProfile_;
		AutocorrelationStats<double> potTempTopProfile_;
		AutocorrelationStats<double> potTempBotProfile_;
		AutocorrelationStats<double> bendistProfile_;
		AutocorrelationStats<double> flowProfile_;
	
    Output(Input const& input);
    void setControlVariates(Input& input, Potential& potential, Galerkin* galerkin);
      
    const double& timeStep() const;
    double& timeStep();
    int& verbose();
    const int& verbose() const;
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
    double energy() const;
    double temperature() const;
    const double& energyMidFlow() const;
    double& energyMidFlow();
		const double& energySumFlow() const;
    double& energySumFlow();
    
    bool doComputeCorrelations() const;
		
		const int& nbOfAutocoPts() const;
		double autocoPtsPeriod() const;
		
    
    ControlVariate& velocityCV();
    ControlVariate& forceCV();
    ControlVariate& lengthCV();
    ControlVariate& midFlowCV();
		ControlVariate& sumFlowCV();
    
    /*Vector<double>& responseForces();
    Vector<double> const& responseForces() const;
    double& responseForces(const int& i);
    const double& responseForces(const int& i) const;*/
    const size_t& decorrelationNbOfIterations() const;
    size_t& decorrelationNbOfIterations();
    double decorrelationTime() const;
		
		void writeProfile(ofstream& out_, size_t iOfIteration);
    //double& integratedAutocorrelationP();
    void display(vector<Particle> const& configuration, size_t iOfIteration);
		void displayGeneratorOnBasis(ofstream& out, vector<Particle> const& configuration, ControlVariate& controlVariate, double time);
    void finalDisplayAutocorrelations();
    void finalDisplay(vector<Particle> const& configuration, Vector<double> const& externalForce);
    void displayFinalFlow(double temperature, double delta_temperature, double tau = nan(""), double xi = 0);
    void displayFinalVelocity(double temperature, double externalForce, int nbOfFourier = 0, int nbOfHermite = 0);

    
    void updateControlVariate(vector<Particle> const& configuration);
		void appendKinTempProfile(double value, size_t iOfIteration, size_t iOfParticle);
		void appendPotTempTopProfile(double value, size_t iOfIteration, size_t iOfParticle);
		void appendPotTempBotProfile(double value, size_t iOfIteration, size_t iOfParticle);
		void appendBendistProfile(double value, size_t iOfIteration, size_t iOfParticle);
		void appendFlowProfile(double value, size_t iOfIteration, size_t iOfParticle);
	};

}
#endif
