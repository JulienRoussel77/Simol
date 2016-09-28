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
  
    //------------- Langevin --------------------

  void writeFinalOutput(Langevin const& /*dyna*/, Isolated const& syst, Output& output)
  {
    if (output.doComputeCorrelations())
      output.finalDisplayCorrelations();
    output.displayFinalVelocity(0, syst.externalForce(0), output.obsVelocity().nbOfFourier(), output.obsVelocity().nbOfHermite());
  }




  //------------- DPDE --------------------

  void samplePositions(DPDE& dyna, Isolated& syst)
  {
    cout << " - Sampling the positions..." << endl;
    syst.getParticle(0).position(0) = 0;
    syst.getParticle(0).momentum(0) = -1.41787;
    syst.getParticle(0).internalEnergy() = 0.0696405;  // TO DO : sample according to equilibrium law?
    dyna.rejectionCount() = 0;
  }



}
