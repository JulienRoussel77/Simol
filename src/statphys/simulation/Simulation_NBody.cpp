#include "Simulation.hpp"

void sampleSystem(Dynamics& dyna, NBody& syst)
{
  //-- initialization of the momenta according to a Gaussian distribution --
  //syst.getParticle(0).momentum() = syst.drawMomentum(dyna.beta(), syst.getParticle(0).mass());
  //-- initialization on a cubic lattice --
  //syst.getParticle(0).position(0) = 0;
}
  
void simulate(Hamiltonian& dyna, NBody& syst)
{
  
}

void computeOutput(Hamiltonian const& /*dyna*/, NBody const& syst, Output& output, size_t /*iOfIteration*/)
{
  
}

//--- CONFLIT DE TEMPLETAGE ICI AUSSI : entre dynamics et system... on ne peut pas preciser que le systeme ?! ---
void writeOutput(Hamiltonian const& /*dyna*/, NBody const& syst, Output& output, size_t iOfIteration)
{
  if (output.doOutput(iOfIteration))
    output.displayObservablesDPDE(syst.configuration(), iOfIteration);
  if (output.doProfileOutput(iOfIteration))
    output.displayParticlesXMakeMol(syst.configuration(), iOfIteration);        
}

void writeFinalOutput(Dynamics const& dyna, NBody const& syst, Output& output)
{
  //output.finalDisplay(syst.configuration(), dyna.externalForce());
}
