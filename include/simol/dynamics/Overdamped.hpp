#ifndef SIMOL_OVERDAMPED_HPP
#define SIMOL_OVERDAMPED_HPP

#include "Dynamics.hpp"

namespace simol
{
  class Overdamped : public Dynamics
  {
    public:
      Overdamped(Input const& input);
      virtual string dynamicsName() const {return "Overdamped";}
      virtual void simulate (System& syst) const;
      //virtual void updateBefore(Particle& particle);
      virtual void updatePosition(Particle& particle) const;
      virtual void getThermo(Output& output) const;
      virtual void getPressure(Output& output) const;
  };

}


#endif
