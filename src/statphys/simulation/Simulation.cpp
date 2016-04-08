#include "Simulation.hpp"

namespace simol {
 
  void sampleSystem(Dynamics& dyna, System& syst)
  {
    dyna.printName();
    syst.printName();
    throw std::invalid_argument("sampleSystem : Function undefined");
  }
  
  void simulate(Dynamics& dyna, System& syst)
  {
    for (auto&& particle : syst.configuration())
      dyna.updateBefore(particle);
    syst.computeAllForces();
    for (auto&& particle : syst.configuration())
      dyna.updateAfter(particle);
  }

  //------------------- writeOutput and its specifications by dynamics ---------------------------

  void writeOutput(Dynamics const& /*dyna*/, System const& /*syst*/, Output& /*output*/, int /*iOfIteration*/)
  {
    throw std::invalid_argument("writeOutput not defined in the general case");
  }

  void writeOutput(Hamiltonian const& /*dyna*/, System const& syst, Output& output, int iOfIteration)
  {
    if (output.doOutput(iOfIteration))
      output.displayObservables(iOfIteration);
    
    if (output.doProfileOutput(iOfIteration))
      output.displayParticles(syst.configuration(), iOfIteration);
  }
  
  void writeOutput(Langevin const& /*dyna*/, System const& syst, Output& output, int iOfIteration)
  {
    if (output.doOutput(iOfIteration))
      output.displayObservables(iOfIteration);
    
    if (output.doProfileOutput(iOfIteration))
      output.displayParticles(syst.configuration(), iOfIteration);
  }

  void writeOutput(DPDE const& /*dyna*/, System const& syst, Output& output, int iOfIteration)
  {
    if (output.doOutput(iOfIteration))
      output.displayObservablesDPDE(syst.configuration(), iOfIteration);
  }

  //---------------------------------------------------------------

  void writeFinalOutput(Dynamics const& /*dyna*/, System const& syst, Output& output)
  {
    output.finalDisplay(syst.configuration(), syst.externalForce());
    if (output.doComputeCorrelations())
      output.finalDisplayAutocorrelations();
  }
   
  
}
