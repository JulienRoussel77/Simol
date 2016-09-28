#include "simol/statphys/simulation/Simulation.hpp"

namespace simol
{
  //void samplePositions(BoundaryLangevin& /*dyna*/, Chain& /*syst*/)
  //{
  //  throw std::invalid_argument("samplePositions(BoundaryLangevin& dyna, Chain& syst) not defined");
  //}

  void writeOutput(BoundaryLangevin const& /*dyna*/, Chain const& syst, Output& output, long int iOfStep)
  {
    if (output.doOutput(iOfStep))// && iOfStep >= 100)
    {
      output.displayThermoVariables(iOfStep);
      output.displayChainPositions(syst.configuration(), iOfStep);
      output.displayChainMomenta(syst.configuration(), iOfStep);

      output.obsMidFlow().display(iOfStep);
      output.obsSumFlow().display(iOfStep);
      output.obsModiFlow().display(iOfStep);
    }

    if (output.doLongPeriodOutput(iOfStep))
    {
      output.displayProfile(iOfStep);
      output.displayParticles(syst.configuration(), iOfStep);
    }
  }


  //-------------- BiChain -----------------

  void samplePositions(BoundaryLangevin& dyna, BiChain& syst)
  {
    cout << " - Sampling the positions..." << endl;
    double alpha, localTemp, localDist;
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
    {
      alpha = iOfParticle / (double) syst.nbOfParticles();
      localTemp = (1 - alpha) * dyna.temperatureLeft() + alpha * dyna.temperatureRight();
      localDist = syst.drawPotLaw(1 / localTemp);
      double prevPosition = (iOfParticle > 0) ? syst.getParticle(iOfParticle - 1).position(0) : 0;

      syst.getParticle(iOfParticle).position(0) = prevPosition + localDist;
    }
  }

  void writeFinalOutput(BoundaryLangevin const& dyna, BiChain const& syst, Output& output)
  {
    output.finalChainDisplay(syst.configuration(), syst.externalForce());
    if (output.doComputeCorrelations())
      output.finalDisplayCorrelations();
    output.displayFinalFlow(dyna.temperature(), dyna.deltaTemperature(), syst.potParameter1(), syst.potParameter2());
  }

  // ------------ TriChain -----------------


  void samplePositions(BoundaryLangevin& dyna, TriChain& syst)
  {
    cout << " - Sampling the positions..." << endl;
    if (syst.isOfFixedVolum())
    {
      cout << "    - Simulation of fixed volum : q_i = 0..." << endl;
      for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
        syst.getParticle(iOfParticle).position(0) = 0;
    }
    else
    {
      double alpha, localTemp, localBending;
      for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      {
        alpha = iOfParticle / (double) syst.nbOfParticles();
        localTemp = (1 - alpha) * dyna.temperatureLeft() + alpha * dyna.temperatureRight();

        localBending = syst.drawPotLaw(1 / localTemp);
        double position1 = (iOfParticle > 0) ? syst.getParticle(iOfParticle - 1).position(0) : 0;
        double position2 = (iOfParticle > 1) ? syst.getParticle(iOfParticle - 2).position(0) : 0;
        syst.getParticle(iOfParticle).position(0) = -position2 + 2 * position1 + localBending;
        syst.getParticle(iOfParticle).momentum() = syst.drawMomentum(1 / localTemp, syst.getParticle(iOfParticle).mass());
      }
    }
  }
  
  void thermalize(Dynamics& /*dyna*/, Chain& /*syst*/) {}

  //void Chain::computeAllForces(Dynamics const& /*model*/)
  //{throw std::invalid_argument("computeAllForces : Function undefined");};*/

  void thermalize(LangevinBase& dyna, Chain& syst)
  {
    //for (auto&& particle : configuration_)
    //dyna.updateBefore(particle);
    for (int i = 0; i < syst.nbOfParticles(); i++)
      dyna.verletFirstPart(syst.getParticle(i));

    syst.computeAllForces();

    for (int i = 0; i < syst.nbOfParticles(); i++)
      dyna.verletSecondPart(syst.getParticle(i));

    for (int i = 0; i < syst.nbOfParticles(); i++)
    {
      double localTemperature = dyna.temperatureLeft() + i * dyna.deltaTemperature() / syst.nbOfParticles();
      dyna.updateOrsteinUhlenbeck(syst.getParticle(0), 1 / localTemperature);
    }
  }
  
  
  
  
  void writeFinalOutput(BoundaryLangevin const& dyna, TriChain const& syst, Output& output)
  {
    output.finalChainDisplay(syst.configuration(), syst.externalForce());
    if (output.doComputeCorrelations())
      output.finalDisplayCorrelations();
    output.displayFinalFlow(dyna.temperature(), dyna.deltaTemperature(), dyna.tauBending(), dyna.xi());
  }



}
