#include "simol/statphys/simulation/Simulation.hpp"

namespace simol
{
  System* createSystem(Input const& input)
  {
    if (input.systemName() == "Isolated")
      return new Isolated(input);
    else if (input.systemName() == "BiChain")
      return new BiChain(input);
    else if (input.systemName() == "TriChain")
      return new TriChain(input);
    else if (input.systemName() == "NBody")
      return new NBody(input);
    else if (input.systemName() == "Colloid")
      return new Colloid(input, 2);
    else if (input.systemName() == "Bicolor")
      return new Bicolor(input);
    else 
      throw std::runtime_error(input.systemName() + " is not a valid system name !");
  }
  
  /*void thermalize(Dynamics& model, System& syst)
  {
    simulate(model, syst);
  }*/
   
  
  


}
