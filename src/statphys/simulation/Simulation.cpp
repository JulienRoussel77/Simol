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
  
  //------------- Overdamped --------------------
  void sampleMomenta(Overdamped& /*dyna*/, System& /*syst*/) {}

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
  
  void sampleMomenta(BoundaryLangevin& dyna, System& syst)
  {
    cout << " - Sampling the momenta..." << endl;
    double alpha, localTemp;
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
    {
      alpha = iOfParticle / (double) syst.nbOfParticles();
      localTemp = (1 - alpha) * dyna.temperatureLeft() + alpha * dyna.temperatureRight();
      syst.getParticle(iOfParticle).momentum() = syst.drawMomentum(1 / localTemp, syst.getParticle(iOfParticle).mass());

      dyna.initializeCountdown(syst.getParticle(iOfParticle));
    }
  }
  
  void thermalize(Dynamics& /*model*/, System& /*syst*/)
  {throw std::invalid_argument("thermalize not defined");}
  
  void thermalize(Dynamics& /*model*/, Isolated& syst)
  {}
  
  
  
  
  void simulate(Dynamics& dyna, System& syst){}
  
  //------------- Hamiltonian -------------------
  void simulate(Hamiltonian& dyna, System& syst)
  {
    for (auto && particle : syst.configuration())
      dyna.verletFirstPart(particle);
    syst.computeAllForces();
    for (auto && particle : syst.configuration())
      dyna.verletSecondPart(particle);
  }
  
  //------------- Overdamped --------------------
  
  void simulate(Overdamped& dyna, System& syst)
  {
    syst.computeAllForces();
    for (auto && particle : syst.configuration())
      dyna.updatePosition(particle);
  }
  
  //------------- LangevinBase --------------------
  
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
  
  //------------ BoundaryLangevin ---------------
  
  void simulate(BoundaryLangevin& dyna, System& syst)
  {
    for (auto && particle: syst.configuration())
      dyna.verletFirstPart(particle);

    syst.computeAllForces();

    for (auto && particle: syst.configuration())
      dyna.verletSecondPart(particle);

    dyna.updateOrsteinUhlenbeck(syst.getParticle(0), dyna.betaLeft());
    dyna.updateOrsteinUhlenbeck(syst.getParticle(syst.nbOfParticles() - 1), dyna.betaRight());

    if (dyna.doMomentaExchange())
      for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles() - 1; iOfParticle++)
        dyna.updateMomentaExchange(syst.getParticle(iOfParticle), syst.getParticle(iOfParticle + 1));
  }
  
  //------------- DPDE --------------------
  
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
  
  
  
  void computeOutput(Dynamics const& dyna, System const& syst, Output& output, long int iOfStep)
  {
    
    // #### faire un output.obsMidFlow().updateBoundaryLangevin(syst.configuration(), dyna.betaLeft(), dyna.betaRight(), dyna.gamma()); ####
    if (output.obsKineticEnergy_) syst.computeKineticEnergy(output);
    if (output.obsPotentialEnergy_) syst.computePotentialEnergy(output);
    syst.computeProfile(output, dyna, iOfStep);
    if (output.obsPressure_) syst.computePressure(output, dyna);
    if (output.obsInternalEnergy_) syst.computeInternalEnergy(output);
    if (output.obsInternalTemperature_) syst.computeInternalTemperature(output, dyna);
    dyna.getThermo(output);
    
    if (output.hasControlVariate()) computeControlVariate(dyna, syst.configuration(), output);
    
    for (auto&& observable : output.observables())
      observable->appendCurrent(iOfStep);
    
    //compute the rejection rate and the negative energies in DPDE
    dyna.specificComputeOutput(output);
  }
  
  void computeControlVariate(Dynamics const& dyna, vector<Particle> const& configuration, Output& output)
  {
    output.cvBasis_.computeValueBasis(configuration);
    dyna.computeGeneratorOnBasis(output.cvBasis_, configuration);
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
  
  void writeOutput(Overdamped const& /*dyna*/, System const& syst, Output& output, long int iOfStep)
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
