#include "simol/statphys/controlVariate/AllGalerkins.hpp"

namespace simol
{

  Galerkin* createGalerkin(Input const& input)
  {
    if (input.doGalerkinCV())
    {
        if (input.dynamicsName() == "Overdamped")
        {
          if (input.potentialName() == "Sinusoidal") return new PeriodicOverdampedGalerkin(input);
          else if (input.potentialName() == "DoubleWell" || input.potentialName() == "Harmonic") return new ColloidOverdampedGalerkin(input);
          else throw runtime_error("This potential matches no Galerkin method !");
        }
        else if (input.dynamicsName() == "Langevin")
        {
          if (input.potentialName() == "Sinusoidal") return new PeriodicLangevinGalerkin(input);
          else if (input.potentialName() == "DoubleWell" || input.potentialName() == "Harmonic") return new ColloidLangevinGalerkin(input);
          else throw runtime_error("This potential matches no Galerkin method !");
        }
        else if (input.dynamicsName() == "BoundaryLangevin") return new BoundaryLangevinGalerkin(input);
        else throw runtime_error("This dynamics matches no Galerkin method !");
      }
    else
      return nullptr;
  }
  
}
