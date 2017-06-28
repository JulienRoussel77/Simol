#ifndef SIMOL_COLLOID_HPP
#define SIMOL_COLLOID_HPP

#include "simol/statphys/system/NBody.hpp"
//#include "simol/statphys/dynamics/LangevinBase.hpp"

namespace simol
{

  class Colloid : public NBody
  {
  public:
    int nbOfColloidParticles_;
    Colloid(Input const& input, int nbOfColloidParticles0);
    int const& nbOfColloidParticles() const;
    int& nbOfColloidParticles();
    virtual void interaction(Particle& particle1, Particle& particle2) const;
    virtual void nonPerInteraction(Particle& particle1, Particle& particle2) const;
    virtual void computeAllForces();
    virtual void samplePositions(DynamicsParameters const& dynaPara);
    
    virtual double length() const;
    virtual double force() const;
  };
  
}

#endif
