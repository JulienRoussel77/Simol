#include "simol/statphys/controlVariate/AllGalerkins.hpp"

namespace simol
{

  Galerkin* createGalerkin(Input const& input, shared_ptr<CVBasis> cvBasis0)
  {
    if (input.doGalerkinCV())
    {
      string potName = input.galerkinPotentialName();
        if (input.dynamicsName() == "Overdamped")
        {
          if (potName == "Sinusoidal") return new PeriodicOverdampedGalerkin(input, cvBasis0);
          else if (potName == "DoubleWell" || potName == "DoubleWellFE" || potName == "Harmonic" || potName == "HarmonicFE" || potName == "TwoTypes") return new ColloidOverdampedGalerkin(input, cvBasis0);
          else throw runtime_error("This potential matches no Galerkin method !");
        }
        else if (input.dynamicsName() == "Langevin")
        {
          if (potName == "Sinusoidal") return new PeriodicLangevinGalerkin(input, cvBasis0);
          else if (potName == "DoubleWell" || potName == "DoubleWellFE" || potName == "Harmonic" || potName == "HarmonicFE") return new ColloidLangevinGalerkin(input, cvBasis0);
          else throw runtime_error("This potential matches no Galerkin method !");
        }
        //else if (input.dynamicsName() == "BoundaryLangevin") return new BoundaryLangevinGalerkin(input);
        else throw runtime_error("This dynamics matches no Galerkin method !");
      }
    else
      return nullptr;
  }
  
}
