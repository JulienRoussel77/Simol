#ifndef SIMOL_OUTPUT_HPP
#define SIMOL_OUTPUT_HPP


#include "tools.hpp"
#include "particle.hpp"
#include "statistics.hpp"
#include "controlVariate.hpp"

namespace simol
{
  


  
  class Output
  {
    std::string outputFoldername_;
    std::ofstream outObservables_;
    std::ofstream outParticles_;
    std::ofstream outReplica_;
		std::ofstream outFinalFlow_;
    std::ofstream outCorrelation_;
    std::ofstream outVelocities_;
    std::ofstream outBeam_;
    std::ofstream outVelocitiesCV_;
    std::ofstream outForcesCV_;
    std::ofstream outLengthsCV_;
    std::ofstream outMidFlowCV_;    
    std::ofstream outMidFlowPT_;
		std::ofstream outSumFlowCV_;    
    std::ofstream outSumFlowPT_;
		std::ofstream outProfile_;
    
    int verbose_;
    size_t periodNumberOfIterations_;
    double timeStep_;
    int dimension_;
    size_t numberOfParticles_;
		size_t numberOfIterations_;
    
    double kineticEnergy_;
    double potentialEnergy_;
    double energyMidFlow_;
		double energySumFlow_;
    
    size_t decorrelationNumberOfIterations_;
    
  public:
    
    ControlVariate* velocityCV_;
    ControlVariate* forceCV_;
    ControlVariate* lengthCV_;
    ControlVariate* midFlowCV_;
		ControlVariate* sumFlowCV_;
		
		AutocorrelationStats<double> temperatureProfile_;
		AutocorrelationStats<double> bendingProfile_;
		AutocorrelationStats<double> flowProfile_;
	
    Output(Input const& input);
    
    void reset(Input const& input, Potential& potential, size_t indexOfReplica);
      
    const double& timeStep() const;
    double& timeStep();
    int& verbose();
    const int& verbose() const;
    double period() const;
    const size_t& periodNumberOfIterations() const;
    bool doOutput(size_t indexOfIteration) const;
		const size_t& numberOfParticles() const;
		const size_t& numberOfIterations() const;
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
    
    ControlVariate* velocityCV();
    ControlVariate* forceCV();
    ControlVariate* lengthCV();
    ControlVariate* midFlowCV();
		ControlVariate* sumFlowCV();
    
    const size_t& decorrelationNumberOfIterations() const;
    size_t& decorrelationNumberOfIterations();
    double decorrelationTime() const;
    //double& integratedAutocorrelationP();
    void display(vector<Particle> const& configuration, size_t indexOfIteration);
    void finalDisplayAutocorrelations();
    void finalDisplay(vector<Particle> const& configuration, Vector<double> const& externalForce);
    void displayFinalFlow(double temperature, double delta_temperature, double tau);

    
    void updateControlVariate(vector<Particle> const& configuration);
		void appendTemperatureProfile(double value, size_t iOfIteration, size_t iOfParticle);
		void appendBendingProfile(double value, size_t iOfIteration, size_t iOfParticle);
		void appendFlowProfile(double value, size_t iOfIteration, size_t iOfParticle);
	};

}
#endif
