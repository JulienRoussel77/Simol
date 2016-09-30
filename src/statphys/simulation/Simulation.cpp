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
      syst(iOfParticle).momentum() = syst.drawMomentum(dyna.beta(), syst(iOfParticle).mass());
  }
  
  //------------- Overdamped --------------------
  void sampleMomenta(Overdamped& /*dyna*/, System& /*syst*/) {}

  //-- initialization of the momenta according to a Gaussian distribution --
  void sampleMomenta(LangevinBase& dyna, System& syst)
  {
    cout << " - Sampling the momenta..." << endl;
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
    {
      syst(iOfParticle).momentum() = syst.drawMomentum(dyna.beta(), syst(iOfParticle).mass());

      dyna.initializeCountdown(syst(iOfParticle));
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
      syst(iOfParticle).momentum() = syst.drawMomentum(1 / localTemp, syst(iOfParticle).mass());

      dyna.initializeCountdown(syst(iOfParticle));
    }
  }
  
  void thermalize(Dynamics& /*model*/, System& /*syst*/)
  {throw std::invalid_argument("thermalize not defined");}
  
  void thermalize(Dynamics& /*model*/, Isolated& /*syst*/)
  {}
  
  
  
  
  void simulate(Dynamics&, System&){}
  
  //------------- Hamiltonian -------------------
  void simulate(Hamiltonian& dyna, System& syst)
  {
    //for (ParticleIterator it = syst.begin(); it != syst.end(); ++it)
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      dyna.verletFirstPart(syst(iOfParticle));
    syst.computeAllForces();

    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      dyna.verletSecondPart(syst(iOfParticle));
  }
  
  //------------- Overdamped --------------------
  
  void simulate(Overdamped& dyna, System& syst)
  {
    syst.computeAllForces();
    
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      dyna.updatePosition(syst(iOfParticle));
  }
  
  //------------- LangevinBase --------------------
  
  void simulate(LangevinBase& dyna, System& syst)
  {
  //for (ParticleIterator it = syst.begin(); !syst.finished(it); ++it)
  for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      dyna.verletFirstPart(syst(iOfParticle));
    syst.computeAllForces();
  //for (ParticleIterator it = syst.begin(); !syst.finished(it); ++it)
  for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      dyna.verletSecondPart(syst(iOfParticle));
  }
  
  //------------- Langevin --------------------
  
  void simulate(Langevin& dyna, System& syst)
  {

    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      dyna.verletFirstPart(syst(iOfParticle));
    syst.computeAllForces();

    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      dyna.verletSecondPart(syst(iOfParticle));

    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      dyna.updateOrsteinUhlenbeck(syst(iOfParticle), dyna.beta());
  }
  
  //------------ BoundaryLangevin ---------------
  
  void simulate(BoundaryLangevin& dyna, System& syst)
  {
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      dyna.verletFirstPart(syst(iOfParticle));

    syst.computeAllForces();

    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      dyna.verletSecondPart(syst(iOfParticle));

    dyna.updateOrsteinUhlenbeck(syst.getParticle(0), dyna.betaLeft());
    dyna.updateOrsteinUhlenbeck(syst.getParticle(syst.nbOfParticles() - 1), dyna.betaRight());

    if (dyna.doMomentaExchange())
      for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles() - 1; iOfParticle++)
        dyna.updateMomentaExchange(syst(iOfParticle), syst(iOfParticle + 1));
  }
  
  //------------- DPDE --------------------
  
  void simulate(DPDE& dyna, System& syst)
  {
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      dyna.verletFirstPart(syst(iOfParticle));
    syst.computeAllForces();
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      dyna.verletSecondPart(syst(iOfParticle));
    //-- fluctuation/dissipation --
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      //dyna.energyReinjection(syst(iOfParticle));  // integration of p at fixed gamma + energy reinjection
      dyna.metropolizedEnergyReinjection(syst(iOfParticle));  // Metropolis correction of effective dynamics on p 
  }
  
  
  
  void computeOutput(Dynamics const& dyna, System const& syst, Output& output, long int iOfStep)
  {
    // #### faire un output.obsMidFlow().updateBoundaryLangevin(syst.configuration(), dyna.betaLeft(), dyna.betaRight(), dyna.gamma()); ####
    if (output.obsKineticEnergy_) dyna.computeKineticEnergy(output, syst);
    if (output.obsPotentialEnergy_) dyna.computePotentialEnergy(output, syst);
    syst.computeProfile(output, iOfStep);
    if (output.obsPressure_) dyna.computePressure(output, syst);
    if (output.obsInternalEnergy_) dyna.computeInternalEnergy(output, syst);
    if (output.obsInternalTemperature_) dyna.computeInternalTemperature(output, syst);
    dyna.getThermo(output);
    
    if (output.hasControlVariate()) computeControlVariate(dyna, syst.configuration(), output);
    
    //cout << "sumFlow = " << output.obsSumFlow().currentValue() << endl;
    
    for (auto&& observable : output.observables())
      observable->appendCurrent(iOfStep);
    
    //compute the rejection rate and the negative energies in DPDE
    dyna.specificComputeOutput(output);
  }
  
  void computeControlVariate(Dynamics const& dyna, vector<Particle*> const& configuration, Output& output)
  {
    output.cvBasis_.computeValueBasis(configuration);
    dyna.computeGeneratorOnBasis(output.cvBasis_, configuration);
  }
  
  //------------------- writeOutput and its specifications by dynamics ---------------------------

  void writeOutput(Dynamics const&, System const& syst, Output& output, long int iOfStep)
  {
    //throw std::invalid_argument("writeOutput not defined in the general case");
    if (output.doOutput(iOfStep))
    {
      for (int iOfObservable=0; iOfObservable < output.nbOfObservables(); iOfObservable++)
        output.observables(iOfObservable)->display(iOfStep);
      output.displayThermoVariables(iOfStep);
      if (output.doOutChain())
      {
        output.displayChainPositions(syst.configuration(), iOfStep);
        output.displayChainMomenta(syst.configuration(), iOfStep);
      }
    }
    
    if (output.doLongPeriodOutput(iOfStep))
    {
      if (output.doOutParticles()) output.displayParticles(syst.configuration(), iOfStep);
      if (output.doOutXMakeMol()) output.displayXMakeMol(syst.configuration(), iOfStep);   
      if (output.doOutBackUp()) output.displayBackUp(syst.configuration(), iOfStep);
      if (output.doOutChain()) output.displayProfile(iOfStep);
    }
  }


  //--------------------------- final output ------------------------------------

  void writeFinalOutput(System const& syst, Output& output)
  {
    output.finalDisplayCorrelations();    
    if (output.doOutChain()) output.finalChainDisplay();
    if (output.doFinalFlow()) output.displayFinalFlow(syst.potParameter1(), syst.potParameter2());
    if (output.doFinalVelocity()) output.displayFinalVelocity();
  }

}
