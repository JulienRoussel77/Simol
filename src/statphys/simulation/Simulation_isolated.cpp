#include "Simulation.hpp"

namespace simol {

  //------------- Hamiltonian --------------------
  
  void sampleSystem(Hamiltonian& /*dyna*/, Isolated& /*syst*/)
  {}
  
  void writeFinalOutput(Hamiltonian const& dyna, Isolated const& syst, Output& output)
  {
    output.finalDisplay(syst.configuration(), syst.externalForce());
    if (output.doComputeCorrelations())
      output.finalDisplayAutocorrelations();
    output.displayFinalVelocity(0, syst.externalForce(0), output.velocityCV_->nbOfFourier(), output.velocityCV_->nbOfHermite());
  }
 
  //------------- Langevin --------------------

  void sampleSystem(LangevinBase& dyna, Isolated& syst)
  {
    syst.getParticle(0).momentum() = syst.drawMomentum(dyna.beta(), syst.getParticle(0).mass());
    syst.getParticle(0).position(0) = syst.drawPotLaw(dyna.beta());
  }
  
  //------------- Overdamped --------------------
  
  void sampleSystem(Overdamped& dyna, Isolated& syst)
  {
    syst.getParticle(0).position(0) = syst.drawPotLaw(dyna.beta());
  }
   
  //------------- DPDE --------------------

  template <>
  void computeOutput(DPDE const& /*dyna*/, Isolated const& syst, Output& output, int /*iOfIteration*/)
  {
    output.kineticEnergy() = syst.getParticle(0).kineticEnergy();
    output.potentialEnergy() = syst.getParticle(0).potentialEnergy();
    output.internalEnergy() = syst.getParticle(0).internalEnergy();
  }

   void sampleSystem(DPDE& dyna, Isolated& syst)
  {
    syst.getParticle(0).momentum() = syst.drawMomentum(dyna.beta(), syst.getParticle(0).mass());
    syst.getParticle(0).position(0) = 0;
    syst.getParticle(0).internalEnergy() = 1;  // TO DO : sample according to equilibrium law?
  }
  
  void simulate(DPDE& dyna, System& syst)
  {
    for (auto&& particle : syst.configuration())
      dyna.verletFirstPart(particle);
    syst.computeAllForces();
    for (auto&& particle : syst.configuration())
      dyna.verletSecondPart(particle);
    //-- fluctuation/dissipation --
    for (int i=0; i < syst.nbOfParticles(); i++)  
      dyna.energyReinjection(syst.getParticle(i));  // integration of p at fixed gamma + energy reinjection 
  }
  
}
