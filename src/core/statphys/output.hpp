#ifndef SIMOL_OUTPUT_HPP
#define SIMOL_OUTPUT_HPP


#include <iostream>
#include "particle.hpp"
#include "statistics.hpp"
#include "controlVariate.hpp"
using std::vector;

namespace simol
{
  


  
  class Output
  {
    std::string outputFoldername_;
    std::ofstream outObservables_;
    std::ofstream outParticles_;
    std::ofstream outReplica_;
    std::ofstream outCorrelation_;
    std::ofstream outVelocities_;
    std::ofstream outBeam_;
    std::ofstream outVelocitiesCV_;
    std::ofstream outForcesCV_;
    std::ofstream outLengthsCV_;
    std::ofstream outFlowCV_;
    
    std::ofstream outFlowPT_;
    
    int verbose_;
    size_t periodNumberOfIterations_;
    double timeStep_;
    int dimension_;
    size_t numberOfParticles_;
    
    double kineticEnergy_;
    double potentialEnergy_;
    double energyFlow_;
    
    //dvec responseForces_;
    
    //dvec velocityRef_;
    //dvec forceRef_;
    //int indexOfIterationRef_;
    
    size_t decorrelationNumberOfIterations_;
    //double integratedAutocorrelationP_;
    
  public:
    //Statistics autocorrelationV;
    //Statistics autocorrelationF;
    
    ControlVariate* velocityCV_;
    ControlVariate* forceCV_;
    ControlVariate* lengthCV_;
    ControlVariate* flowCV_;
	
    Output(Input const& input);
    
    void reset(Input const& input, Potential& potential, size_t indexOfReplica);
      
    const double& timeStep() const;
    double& timeStep();
    int& verbose();
    const int& verbose() const;
    double period() const;
    const size_t& periodNumberOfIterations() const;
    bool doOutput(size_t indexOfIteration) const;
    
    const double& kineticEnergy() const;
    double& kineticEnergy();
    const double& potentialEnergy() const;
    double& potentialEnergy();
    double energy() const;
    double temperature() const;
    const double& energyFlow() const;
    double& energyFlow();
    
    bool doComputeCorrelations() const;
    
    ControlVariate* velocityCV();
    ControlVariate* forceCV();
    ControlVariate* lengthCV();
    ControlVariate* flowCV();
    
    /*dvec& responseForces();
    dvec const& responseForces() const;
    double& responseForces(const int& i);
    const double& responseForces(const int& i) const;*/
    const size_t& decorrelationNumberOfIterations() const;
    size_t& decorrelationNumberOfIterations();
    double decorrelationTime() const;
    //double& integratedAutocorrelationP();
    void display(vector<Particle> const& configuration, size_t indexOfIteration);
    void finalDisplayAutocorrelations();
    void finalDisplay(vector<Particle> const& configuration, dvec const& externalForce, size_t indexOfIteration);
    //void displayVelocity(size_t indexOfIteration);

    
    void updateControlVariate(vector<Particle> const& configuration);
  };

}
#endif