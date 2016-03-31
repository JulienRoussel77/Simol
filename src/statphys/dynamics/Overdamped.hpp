#ifndef SIMOL_OVERDAMPED_HPP
#define SIMOL_OVERDAMPED_HPP

#include "UniformStochasticDynamics.hpp"

namespace simol
{
  class Overdamped : public UniformStochasticDynamics
  {
    public:
        Overdamped(Input const& input);
        virtual void updateBefore(Particle& particle);
        virtual void updateAfter(Particle& particle);
 };

}


#endif
