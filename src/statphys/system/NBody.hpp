#ifndef SIMOL_NBODY_HPP
#define SIMOL_NBODY_HPP
#include "System.hpp"

namespace simol
{

  class NBody : public System
  {
  public:
    NBody(Input const& input);
    void printName() const;
    void computeAllForces();
    int nbOfParticlesPerDimension() const;
    double latticeParameter() const;
    void interaction(Particle& particle1, Particle& particle2) const;
  protected:
    int nbOfParticlesPerDimension_;
    double latticeParameter_;
    double domainSize_;
    //ofstream DEBUG_;
   };
  
}


#endif