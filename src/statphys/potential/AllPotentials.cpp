#include "simol/statphys/potential/AllPotentials.hpp"

namespace simol
{
  
  Potential* createPotential(Input const& input)
  {
    if (input.potentialName() == "Sinusoidal")
      return new Sinusoidal(input);
    else if (input.potentialName() == "SumSinusoidal")
      return new SumSinusoidal(input);
    else if (input.potentialName() == "FracSinusoidal")
      return new FracSinusoidal(input);
    else if (input.potentialName() == "DoubleWell")
      return new DoubleWell(input);
    else if (input.potentialName() == "Harmonic")
      return new Harmonic(input);
    else if (input.potentialName() == "Rotor")
      return new Rotor(input);
    else if (input.potentialName() == "FPU")
      return new FPU(input);
    else if (input.potentialName() == "SpaceSinus")
      return new SpaceSinus(input);
    else if (input.potentialName() == "LennardJones")
      return new LennardJones(input);
    else
      std::cout << input.potentialName() << " is not a valid potential !" << std::endl;
    
    return 0;
  }
  
}
