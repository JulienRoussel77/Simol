#ifndef SIMOL_ALLGALERKINS_HPP
#define SIMOL_ALLGALERKINS_HPP

#include "Galerkin.hpp"
#include "LangevinGalerkin.hpp"
#include "OverdampedGalerkin.hpp"

namespace simol
{

  Galerkin* createGalerkin(Input const& input);

}

#endif