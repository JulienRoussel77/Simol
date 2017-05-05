#ifndef SIMOL_BICOLOR_HPP
#define SIMOL_BICOLOR_HPP

#include "simol/statphys/system/NBody.hpp"
//#include "simol/statphys/dynamics/LangevinBase.hpp"

namespace simol
{

  class Bicolor : public NBody
  {
    double fixedVelocity_;
  public:
    Bicolor(Input const& input);
    double const& fixedVelocity() const {return fixedVelocity_;}
    virtual void computeAllForces();
    
  };
  
}

#endif
