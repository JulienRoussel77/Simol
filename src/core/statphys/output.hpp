#ifndef SIMOL_OUTPUT_HPP
#define SIMOL_OUTPUT_HPP


#include <iostream>
#include "particle.hpp"

namespace simol
{

  class Output
  {
    std::string outputFoldername_;
    std::ofstream outParticles_;
    std::ofstream outReplica_;
    
    dvec sumForces_;
    std::vector<dvec> responseForces_;
  public:
    Output(Input const& input);
    void initialize();
    dvec& sumForces();
    std::vector<dvec>& responseForces();
    const dvec& responseForces(int const& i) const;
    dvec& responseForces(int const& i);
    void display(Particle const& particle, double time);
    void finalDisplay(Particle const& particle, double time);
    void finalDisplayExternalForce(Particle const& particle, dvec const& externalForce, dvec const& responseForce, double time);
    int verbose_;
  };

}
#endif