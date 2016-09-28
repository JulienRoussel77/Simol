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
    output.displayFinalVelocity(0, syst.externalForce(0), output.obsVelocity().nbOfFourier(), output.obsVelocity().nbOfHermite());
    //vector<double> parameters = {0, syst.externalForce(0)};//, output.velocity().nbOfFourier(), output.velocity().nbOfHermite()};
    //output.obsVelocity().finalDisplay(parameters);
  }
  
  //------------- Overdamped --------------------
    
  /*template <>
  void computeOutput(Overdamped const& dyna, Isolated const& syst, Output& output, long int iOfStep)
  {
    output.potentialEnergy() = syst.getParticle(0).potentialEnergy();
    // In the case of the trichain we add the potential of the wall interaction
    output.potentialEnergy() += syst.boundaryPotEnergy();
    updateAllControlVariates(dyna, syst, output, iOfStep);
  }*/
  
    //------------- Langevin --------------------

  void writeFinalOutput(Langevin const& /*dyna*/, Isolated const& syst, Output& output)
  {
    if (output.doComputeCorrelations())
      output.finalDisplayCorrelations();
    output.displayFinalVelocity(0, syst.externalForce(0), output.obsVelocity().nbOfFourier(), output.obsVelocity().nbOfHermite());
    //vector<double> parameters = {0, syst.externalForce(0)};//, output.velocity().nbOfFourier(), output.velocity().nbOfHermite()};
    //output.obsVelocity().finalDisplay(parameters);
  }

  //------------- Overdamped --------------------

  void sampleMomenta(Overdamped& /*dyna*/, Isolated& /*syst*/) {}


  //------------- DPDE --------------------

  /*template <>
  void computeOutput(DPDE const& dyna, Isolated const& syst, Output& output, long int iOfStep)
  {
    //-- instantaneous values --
    output.kineticEnergy() = syst.getParticle(0).kineticEnergy();
    output.potentialEnergy() = syst.getParticle(0).potentialEnergy();
    output.internalEnergy() = syst.getParticle(0).internalEnergy();
    // -- averages of observables --
    output.appendKineticEnergy(syst.getParticle(0).kineticEnergy(), iOfStep);
    output.appendPotentialEnergy(syst.getParticle(0).potentialEnergy(), iOfStep);
    output.appendInternalEnergy(syst.getParticle(0).internalEnergy(), iOfStep);
    //-- rejection rate --
    output.rejectionCount() = dyna.rejectionCount();
  }*/

  void samplePositions(DPDE& dyna, Isolated& syst)
  {
    cout << " - Sampling the positions..." << endl;
    syst.getParticle(0).position(0) = 0;
    syst.getParticle(0).momentum(0) = -1.41787;
    syst.getParticle(0).internalEnergy() = 0.0696405;  // TO DO : sample according to equilibrium law?
    dyna.rejectionCount() = 0;
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
      //dyna.energyReinjection(syst.getParticle(iOfParticle));  // integration of p at fixed gamma + energy reinjection
      dyna.metropolizedEnergyReinjection(syst.getParticle(iOfParticle));  // Metropolis correction of effective dynamics on p 
  }

}
