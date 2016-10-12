#include "simol/statphys/simulation/Simulation.hpp"

namespace simol
{
  System* createSystem(Input const& input)
  {
    if (input.systemName() == "Isolated")
      return new Isolated(input);
    else if (input.systemName() == "BiChain")
      return new BiChain(input);
    else if (input.systemName() == "TriChain")
      return new TriChain(input);
    else if (input.systemName() == "NBody")
      return new NBody(input);
    else 
      throw std::runtime_error(input.systemName() + " is not a valid system name !");
  }
  

  
  void sampleInternalEnergies(Dynamics const&, System&)
  {}
  
  void sampleInternalEnergies(DPDE const&, System& syst)
  {
    cout << " - Sampling internal energies..." << endl;
    for (int i = 0; i < syst.nbOfParticles(); i++)
      syst.getParticle(i).internalEnergy() = 1;
  }
  
  void thermalize(Dynamics& model, System& syst)
  {
    simulate(model, syst);
  }
   
  
  
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
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      dyna.updatePosition(syst(iOfParticle));
    
    syst.computeAllForces();
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
  
  //------------- DPDE [1D systems...] --------------------
  
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
    // #### faire un output.obsMidFlow().updateBoundaryLangevin(syst, dyna.betaLeft(), dyna.betaRight(), dyna.gamma()); ####
    if (output.obsKineticEnergy_) dyna.computeKineticEnergy(output, syst);
    if (output.obsPotentialEnergy_) dyna.computePotentialEnergy(output, syst);
    if (syst.isBiChain()) dyna.computeProfileBiChain(output, syst, iOfStep);
    else if (syst.isTriChain()) dyna.computeProfileBiChain(output, syst, iOfStep);
    if (output.obsPressure_) dyna.computePressure(output, syst);
    if (output.obsInternalEnergy_) dyna.computeInternalEnergy(output, syst);
    if (output.obsInternalTemperature_) dyna.computeInternalTemperature(output, syst);
    if (output.obsVelocity_) output.velocity() = syst(0).velocity(0);
    if (output.obsForce_) output.force() = syst(0).force(0);
    dyna.getThermo(output);
    
    if (output.hasControlVariate()) computeControlVariate(dyna, syst, output);
    
    //cout << "sumFlow = " << output.obsSumFlow().currentValue() << endl;
    
    for (auto&& observable : output.observables())
      observable->appendCurrent(iOfStep);
    
    //compute the rejection rate and the negative energies in DPDE
    dyna.specificComputeOutput(output);
  }
  
  void computeControlVariate(Dynamics const& dyna, System const& syst, Output& output)
  {
    output.cvBasis_->computeValueBasis(syst);
    dyna.computeGeneratorOnBasis(*output.cvBasis_, syst);
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
        output.displayChainPositions(syst, iOfStep);
        output.displayChainMomenta(syst, iOfStep);
      }
    }
    
    if (output.doLongPeriodOutput(iOfStep))
    {
      if (output.doOutParticles()) output.displayParticles(syst, iOfStep);
      if (output.doOutXMakeMol()) output.displayXMakeMol(syst, iOfStep);   
      if (output.doOutBackUp()) output.displayBackUp(syst, iOfStep);
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
