#ifndef SIMOL_BICOLOR_HPP
#define SIMOL_BICOLOR_HPP

#include "simol/statphys/system/NBody.hpp"
//#include "simol/statphys/dynamics/LangevinBase.hpp"

namespace simol
{

  class Bicolor : public NBody
  {
  public:
    Bicolor(Input const& input);
    virtual void computeAllForces();
    virtual double velocity() const;
    virtual double force() const;
    virtual void enforceConstraint(double& lagrangeMultiplier, double fixedVelocity);
  };
  
}

#endif
