#ifndef SIMOL_OUTPUT_HPP
#define SIMOL_OUTPUT_HPP


#include <iostream>
#include "particle.hpp"
using std::vector;

namespace simol
{
  
  class Statistics
  {
    vector<double> values_;
    vector<int> nbValues_;
  public:
    Statistics(int nbIndices);
    //virtual ~Statistics(){};
    void append(int i, double value);
    double& operator()(int i);
    const double& operator()(int i) const;
  };
  
  class Output
  {
    std::string outputFoldername_;
    std::ofstream outObservables_;
    std::ofstream outParticles_;
    std::ofstream outReplica_;
    std::ofstream outCorrelation_;
    
    int verbose_;
    double timeStep_;
    int dimension_;
    int numberOfParticles_;
    
    double kineticEnergy_;
    double potentialEnergy_;
    
    dvec responseForces_;
    
    dvec velocityRef_;
    dvec forceRef_;
    double timeRef_;
    
    size_t decorrelationNumberOfIterations_;
    //double integratedAutocorrelationP_;
  public:
    Statistics autocorrelationV;
    Statistics autocorrelationF;
	
    Output(Input const& input);
    void initialize(dvec const& refVelocity = 0);
    const double& timeStep() const;
    int& verbose();
    const int& verbose() const;
    const double& kineticEnergy() const;
    double& kineticEnergy();
    const double& potentialEnergy() const;
    double& potentialEnergy();
    double energy() const;
    double temperature() const;
    bool doComputeCorrelations() const;
    dvec& responseForces();
    const dvec& responseForces() const;
    double& responseForces(const int& i);
    const double& responseForces(const int& i) const;
    double& timeRef();
    const double& timeRef() const;
    dvec& velocityRef();
    const dvec& velocityRef() const;
    dvec& forceRef();
    const dvec& forceRef() const;
    const size_t& decorrelationNumberOfIterations() const;
    size_t& decorrelationNumberOfIterations();
    //double& integratedAutocorrelationP();
    void display(vector<Particle> const& configuration, double time);
    void finalDisplayAutocorrelations();
    void finalDisplay(vector<Particle> const& configuration, dvec const& externalForce, double time);
    void appendAutocorrelationV(dvec const& velocity, double time);
    void appendAutocorrelationF(dvec const& force, double time);
  };

}
#endif