#include "simol/statphys/controlVariate/AllGalerkins.hpp"

namespace simol
{

  Galerkin* createGalerkin(Input const& input)
  {
    if (input.doGalerkinCV())
    {
        if (input.dynamicsName() == "Overdamped") return new OverdampedGalerkin(input);
        else if (input.dynamicsName() == "PeriodicLangevin") return new PeriodicLangevinGalerkin(input);
        else if (input.dynamicsName() == "ColloidLangevin") return new ColloidLangevinGalerkin(input);
        else if (input.dynamicsName() == "BoundaryLangevin") return new BoundaryLangevinGalerkin(input);
        else throw runtime_error("This dynamics matches no Galerkin method !");
      }
    else
      return nullptr;
  }
  
}
