#ifndef SIMOL_OVERDAMPED_HPP
#define SIMOL_OVERDAMPED_HPP

#include "Dynamics.hpp"

namespace simol
{
  class Overdamped : public Dynamics
  {
    public:
        Overdamped(Input const& input);
        virtual void updateBefore(Particle& particle);
        virtual void updateAfter(Particle& particle);
 };

}


#endif
