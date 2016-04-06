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

namespace simol
{

  Potential* createPotential(Input const& input);

}

#endif