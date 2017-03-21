#ifndef SIMOL_ALLPOTENTIALS_HPP
#define SIMOL_ALLPOTENTIALS_HPP

#include "DoubleWell.hpp"
#include "FPU.hpp"
#include "FracSinusoidal.hpp"
#include "Harmonic.hpp"
#include "LennardJones.hpp"
#include "Rotor.hpp"
#include "Sinusoidal.hpp"
#include "SpaceSinus.hpp"
#include "SumSinusoidal.hpp"
#include "SoftDPD.hpp"

namespace simol
{
  Potential* createPotential(Input const& input, string potName);
  Potential* createPotential(Input const& input);
  
  class TwoTypes : public Potential
  {
    public:
      TwoTypes(Input const& input);
      double operator()(double position, int type) const;
      DVec gradient(double position, int type) const;
    private:
      Potential* firstPot_;
      Potential* secondPot_;
  };

}

#endif
