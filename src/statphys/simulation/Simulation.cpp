#include "simol/statphys/simulation/Simulation.hpp"

namespace simol
{

  void samplePositions(Dynamics& dyna, System& syst)
  {
    dyna.printName();
    syst.printName();
    throw std::invalid_argument("samplePositions : Function undefined");
  }

  //-- initialization of the momenta according to a Gaussian distribution --
  void sampleMomenta(Dynamics& dyna, System& syst)
  {
    cout << " - Sampling the momenta..." << endl;
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      syst.getParticle(iOfParticle).momentum() = syst.drawMomentum(dyna.beta(), syst.getParticle(iOfParticle).mass());
  }

  //-- initialization of the momenta according to a Gaussian distribution --
  void sampleMomenta(LangevinBase& dyna, System& syst)
  {
    cout << " - Sampling the momenta..." << endl;
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
    {
      syst.getParticle(iOfParticle).momentum() = syst.drawMomentum(dyna.beta(), syst.getParticle(iOfParticle).mass());

      dyna.initializeCountdown(syst.getParticle(iOfParticle));
    }
  }


  /*void simulate(Dynamics& dyna, System& syst)
  {
  for (auto && particle : syst.configuration())
      dyna.updateBefore(particle);
    syst.computeAllForces();
  for (auto && particle : syst.configuration())
      dyna.updateAfter(particle);
  }*/
  
  void simulate(LangevinBase& dyna, System& syst)
  {
  for (auto && particle : syst.configuration())
      dyna.verletFirstPart(particle);
    syst.computeAllForces();
  for (auto && particle : syst.configuration())
      dyna.verletSecondPart(particle);
  }
  
  //------------- Langevin --------------------
  
  void simulate(Langevin& dyna, System& syst)
  {
    for (auto && particle : syst.configuration())
      dyna.verletFirstPart(particle);
    syst.computeAllForces();
    for (auto && particle : syst.configuration())
      dyna.verletSecondPart(particle);
    for (auto && particle : syst.configuration())
      dyna.updateOrsteinUhlenbeck(particle, dyna.beta());
  }

  //------------------- writeOutput and its specifications by dynamics ---------------------------

  void writeOutput(Dynamics const& /*dyna*/, System const& /*syst*/, Output& /*output*/, long int /*iOfStep*/)
  {
    throw std::invalid_argument("writeOutput not defined in the general case");
  }

  void writeOutput(Hamiltonian const& /*dyna*/, System const& syst, Output& output, long int iOfStep)
  {
    if (output.doOutput(iOfStep))
      output.displayThermoVariables(iOfStep);

    if (output.doLongPeriodOutput(iOfStep))
      output.displayParticles(syst.configuration(), iOfStep);
  }

  void writeOutput(Langevin const& dyna, System const& syst, Output& output, long int iOfStep)
  {
    if (output.doOutput(iOfStep))
    {
      output.displayThermoVariables(iOfStep);
      output.obsVelocity().display(iOfStep);
    }
    if (output.doLongPeriodOutput(iOfStep))
      output.displayParticles(syst.configuration(), iOfStep);
  }

  void writeOutput(DPDE const& /*dyna*/, System const& syst, Output& output, long int iOfStep)
  {
    if (output.doOutput(iOfStep))
      output.displayThermoVariablesDPDE(syst.configuration(), iOfStep);
  }

  //--------------------------- final output ------------------------------------

  void writeFinalOutput(Dynamics const& /*dyna*/, System const& syst, Output& output)
  {
    if (output.doComputeCorrelations())
      output.finalDisplayCorrelations();
      
  }

  void writeFinalOutput(DPDE const& /*dyna*/, System const& /*syst*/, Output& output)
  {
    if (output.doComputeCorrelations())
      output.finalDisplayCorrelationsDPDE();
  }


}
