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
      //virtual void updateBefore(Particle& particle);
      virtual void updatePosition(Particle& particle);
      virtual void getThermo(Output& output) const;
      virtual void getPressure(Output& output) const;
      virtual void computeGeneratorOnBasis(CVBasis& cvBasis, System const& syst) const;
  };

}


#endif
