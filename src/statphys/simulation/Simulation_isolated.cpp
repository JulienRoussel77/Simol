#include "simol/statphys/simulation/Simulation.hpp"

namespace simol
{

  void samplePositions(Dynamics& dyna, Isolated& syst)
  {
    cout << " - Sampling the positions..." << endl;
    syst.getParticle(0).position(0) = syst.drawPotLaw(dyna.beta());
  }

  //------------- Hamiltonian --------------------

  void writeFinalOutput(Hamiltonian const& /*dyna*/, Isolated const& syst, Output& output)
  {
    if (output.doComputeCorrelations())
      output.finalDisplayCorrelations();
    output.displayFinalVelocity(0, syst.externalForce(0), output.velocityCV_->nbOfFourier(), output.velocityCV_->nbOfHermite());
  }

  //------------- Overdamped --------------------

  void sampleMomenta(Overdamped& /*dyna*/, Isolated& /*syst*/) {}


  //------------- DPDE --------------------

  template <>
  void computeOutput(DPDE const& /*dyna*/, Isolated const& syst, Output& output, long int iOfStep)
  {
    //-- instantaneous values --
    output.kineticEnergy() = syst.getParticle(0).kineticEnergy();
    output.potentialEnergy() = syst.getParticle(0).potentialEnergy();
    output.internalEnergy() = syst.getParticle(0).internalEnergy();
    // -- averages --
    output.appendKineticEnergy(syst.getParticle(0).kineticEnergy(), iOfStep);
  }

  void samplePositions(DPDE& /*dyna*/, Isolated& syst)
  {
    cout << " - Sampling the positions..." << endl;
    syst.getParticle(0).position(0) = 0;
    syst.getParticle(0).internalEnergy() = 1;  // TO DO : sample according to equilibrium law?
  }

  void simulate(DPDE& dyna, System& syst)
  {
    for (auto && particle : syst.configuration())
      dyna.verletFirstPart(particle);
    syst.computeAllForces();
    for (auto && particle : syst.configuration())
      dyna.verletSecondPart(particle);
    //-- fluctuation/dissipation --
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      dyna.energyReinjection(syst.getParticle(iOfParticle));  // integration of p at fixed gamma + energy reinjection
  }

}
