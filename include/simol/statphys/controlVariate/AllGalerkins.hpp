#ifndef SIMOL_ALLGALERKINS_HPP
#define SIMOL_ALLGALERKINS_HPP

#include "Galerkin.hpp"
#include "LangevinGalerkin.hpp"
#include "OverdampedGalerkin.hpp"

namespace simol
{

  Galerkin* createGalerkin(const Input& input, shared_ptr<CVBasis> cvBasis0);

}

#endif
