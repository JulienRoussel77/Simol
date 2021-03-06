//#ifndef SIMOL_DYNAMICS_IPP
//#define SIMOL_DYNAMICS_IPP

#include "simol/dynamics/Dynamics.hpp"

using std::cout;
using std::endl;

namespace simol
{



  ///
  ///Constructor
  Dynamics::Dynamics(Input const& input):
    parameters_(input),
    timeStep_(input.timeStep()),
    nbOfSteps_(input.nbOfSteps()),
    thermalizationNbOfSteps_(input.thermalizationNbOfSteps()),
    burninNbOfSteps_(input.burninNbOfSteps())
  {}




  ///
  ///Read-write accessor for the time step "dt"
  double& Dynamics::timeStep() {return timeStep_;}
  ///
  ///Read-only accessor for the time step "dt"
  double const& Dynamics::timeStep() const {return timeStep_;}
  ///
  ///Read-write accessor for the number of steps of the simulation
  long int& Dynamics::nbOfSteps() {return nbOfSteps_;}
  ///
  ///Read-only accessor for the number of steps of the simulation
  long int const& Dynamics::nbOfSteps() const {return nbOfSteps_;}
  ///
  ///Returns the total time "t_f" of the simulation
  double Dynamics::finalTime() const {return timeStep_ * nbOfSteps_;}
  ///
  ///Read-write accessor for the number of steps of the thermalization
  long int& Dynamics::thermalizationNbOfSteps() {return thermalizationNbOfSteps_;}
  ///
  ///Read-only accessor for the number of steps of the thermalization
  long int const& Dynamics::thermalizationNbOfSteps() const {return thermalizationNbOfSteps_;}
  ///
  ///Read-write accessor for the number of steps of the burnIn
  long int& Dynamics::burninNbOfSteps() {return burninNbOfSteps_;}
  ///
  ///Read-only accessor for the number of steps of the burnIn
  long int const& Dynamics::burninNbOfSteps() const {return burninNbOfSteps_;}

  const std::shared_ptr<RNG>& Dynamics::rng() const {return rng_;}

  std::shared_ptr<RNG>& Dynamics::rng() {return rng_;}

  ///
  ///Returns the temperature
  double const& Dynamics::temperature() const {return parameters_.temperature();}
  ///
  ///Read-only access for the temperature
  double const& Dynamics::temperatureLeft() const {return parameters_.temperatureLeft();}
  ///
  ///Read-only access for the temperature
  double const& Dynamics::temperatureRight() const {return parameters_.temperatureRight();}
  ///Returns the difference between the mean temperature and the one at the left end
  ///This is equal to {eta}
  double Dynamics::deltaTemperature() const
  {
    return (temperatureLeft() - temperatureRight()) / 2;
  }
  ///
  ///Read-only access for the inverse temperature
  double const& Dynamics::beta() const {return parameters_.beta();}
  ///
  ///Read-only access for the inverse temperature
  double const& Dynamics::betaLeft() const {return parameters_.betaLeft();}
  ///
  ///Read-only access for the inverse temperature
  double const& Dynamics::betaRight() const {return parameters_.betaRight();}

  //bool const& Dynamics::isMollified() const {return parameters_.isMollified();}
  

  //----------- Ouputs ------------------
  
  ///
  /// Compute observables such as energies, thermodynamic quantities and linear combinations of velocities
  /// Is followed by writeOutput
  void Dynamics::computeOutput(System const& syst, Output& output, long int iOfStep) const
  {
    if (output.obsKineticEnergy_) computeKineticEnergy(output, syst);
    if (output.obsPotentialEnergy_) computePotentialEnergy(output, syst);
    if (syst.isBiChain()) computeProfileBiChain(output, syst, iOfStep);
    else if (syst.isTriChain()) computeProfileBiChain(output, syst, iOfStep);
    if (output.obsPressure_) computePressure(output, syst);
    if (output.obsInternalEnergy_) computeInternalEnergy(output, syst);
    if (output.obsInternalTemperature_) computeInternalTemperature(output, syst);
    if (output.obsLength_) output.length() = syst.length();
    if (output.obsVelocity_) output.velocity() = syst.velocity();
    if (output.obsForce_) output.force() = syst.force();
    if (output.obsLagrangeMultiplier_) output.lagrangeMultiplier() = syst.lagrangeMultiplier();
    getThermo(output);
    
       
    for (auto&& observable : output.observables())
      observable->appendCurrent(iOfStep);
    
    //compute specific outputs, e.g. rejection rate, negative energies in DPDE, etc
    specificComputeOutput(output);
  }
  
  //------------------- writeOutput and its specifications by dynamics ---------------------------
  
  ///
  /// Write outputs in the output files for relevent observables, some of these computed in computeOutput
  void Dynamics::writeOutput(System const& syst, Output& output, long int iOfStep) const
  {
    // Test if the iteration index is a multiple of the input parameter "PrintPeriod"
    if (output.doOutput(iOfStep))
    {
      for (int iOfObservable=0; iOfObservable < output.nbOfObservables(); iOfObservable++)
        output.observables(iOfObservable)->display(iOfStep);
      output.displayThermoVariables(iOfStep);
      if (output.doOutChain())
        output.displayInstantProfile(syst, iOfStep);
    }
    
    // Test if the iteration index is a multiple of the input parameter "PrintLongPeriod"
    // These outut are more time and memory consuming so they are written less often
    if (output.doLongPeriodOutput(iOfStep))
    {
      if (output.doOutParticles()) output.displayParticles(syst, iOfStep);
      if (output.doXMakeMol()) output.displayXMakeMol(syst, iOfStep);   
      if (output.doOutBackUp()) output.displayBackUp(syst, iOfStep);
      if (output.doOutChain()) output.displayProfile(iOfStep);
    }
  }


  //--------------------------- final output ------------------------------------

  ///
  /// Similar to writeOutput but called at the end of the simulation for one-line summaries
  void Dynamics::writeFinalOutput(System const& syst, Output& output) const
  {
    output.finalDisplayCorrelations();    
    if (output.doOutChain()) output.finalChainDisplay();
    if (output.doFinalFlux()) output.displayFinalFlux(syst.potParameter1(), syst.potParameter2(), syst.pairPotential().harmonicFrequency());
    if (output.doFinalLength()) output.displayFinalLength();
    if (output.doFinalVelocity()) output.displayFinalVelocity();
    if (output.doFinalLagrangeMultiplier())
    {
      if (output.doOutChain())
        output.displayFinalChainLagrangeMultiplier(syst.potParameter1(), syst.potParameter2());
      else
        output.displayFinalLagrangeMultiplier();        
    }
  }
  
  /// The simulation starts with a "ThermalizationTime" where an iteration is described by this function
  /// By default this is a standard iteration but can be replaced by a smarter dynamics, with a shorter transitory regime
  void Dynamics::thermalize(System& syst) const
  {
    simulate(syst);
  }
  
  /// Initializes the system, either using a configuration file prescribing all initial variables or a distribution law
  /// Ideally the inital state is sampled according to the invariant probability measure of the synamics
  void Dynamics::sampleSystem(System& syst) const
  {
    cout << " Initialization of the system..." << endl;

    // Test if a configuration file has been provided, if not we use distribution law
    if (!syst.doSetting())
    {
      syst.sampleMomenta(parameters());
      syst.samplePositions(parameters());
      syst.sampleInternalEnergies();
    }
    
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      syst(iOfParticle).oldGaussian() = syst.rng()->gaussian();

    
    syst.computeAllForces();
    cout << " - Thermalization (" << thermalizationNbOfSteps() << " steps)..." << endl;

    for (long int iOfStep  = 0; iOfStep < thermalizationNbOfSteps(); ++iOfStep)
      thermalize(syst);

    cout << " - Burn-in (" << burninNbOfSteps() << " steps)..." << endl;

    for (long int iOfStep  = 0; iOfStep < burninNbOfSteps(); ++iOfStep)
      simulate(syst);
    
    cout << " Starting production mode" << endl;
    cout << endl;
  }
  
  ///
  /// Main function
  void Dynamics::launch(System& syst, Output& output)
  {    
    //---- initialization (including burn-in) -----
    sampleSystem(syst);
    //---- actual steps -----
    for (long int iOfStep  = 0; iOfStep < nbOfSteps(); ++iOfStep)
    {
      //--- display progress every time 10% of simulation elapsed ---
      if ((10 * iOfStep) % nbOfSteps() == 0)
        cout << "---- Run " << (100 * iOfStep) / nbOfSteps() << " % completed ----" << endl;

      //--- compute and write outputs if required ----
      computeOutput(syst, output, iOfStep);
      writeOutput(syst, output, iOfStep);
      //---- update the system_ by the numerical integration ---
      simulate(syst);
    }

    //--- write final outputs ----
    writeFinalOutput(syst, output);
  }

  ///
  ///Standard first part of the numerical integration : half upadate on "p" and update on "q"
  void Dynamics::verletFirstPart(Particle& particle) const
  {
    particle.momentum() += timeStep_ * particle.force() / 2;
    particle.position() += timeStep_ * particle.momentum() / particle.mass();
  }

  ///
  ///Standard second part of the numerical integration : half upadate on "p"
  void Dynamics::verletSecondPart(Particle& particle) const
  {
    particle.momentum() += timeStep_ * particle.force() / 2;
  }
  
  ///
  /// Integrates dp = f(q) dt over half a timestep for one particle
  void Dynamics::updateMomentum(Particle& particle) const
  {
    particle.momentum() += timeStep_ * particle.force() / 2;
  }
  
  ///
  /// Integrates dq = p/m dt over half a timestep for one particle
  void Dynamics::updatePosition(Particle& particle) const
  {
    particle.position() += timeStep_ * particle.momentum() / particle.mass();
  }
  
  // --- computes thermodynamic quantities and stock the result in observable objects ---
  
  void Dynamics::computeKineticEnergy(Output& output, System const& syst) const
  {
    output.obsKineticEnergy().currentValue() = 0;
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      output.obsKineticEnergy().currentValue() += syst(iOfParticle).kineticEnergy();
  }
  
  void Dynamics::computePotentialEnergy(Output& output, System const& syst) const
  {
    output.obsPotentialEnergy().currentValue() = 0;
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      output.obsPotentialEnergy().currentValue() += syst(iOfParticle).potentialEnergy();
    output.obsPotentialEnergy().currentValue() += syst.boundaryPotEnergy();
  }
  
  void Dynamics::computePressure(Output& output, System const& syst) const
  {
    output.totalVirial() = 0;
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      output.totalVirial() += syst(iOfParticle).virial();
    
    //Computes the instantaneous pressure, knowing the total virial
    getPressure(output);
  }
  
  void Dynamics::computeInternalEnergy(Output& output, System const& syst) const
  {
    output.obsInternalEnergy().currentValue() = 0;
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      output.obsInternalEnergy().currentValue() += syst(iOfParticle).internalEnergy();
  }
  
  void Dynamics::computeInternalTemperature(Output& output, System const& syst) const
  {
    output.obsInternalTemperature().currentValue() = 0;
    //--- version for "general microEOS" where exp(s(eps)) does not vanish at eps = 0 --------
    //for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
    //  output.obsInternalTemperature().currentValue() += syst(iOfParticle).internalEnergy()*entropy_derivative(syst(iOfParticle).internalEnergy()); 
    //output.obsInternalTemperature().currentValue() = output.obsInternalTemperature().currentValue()/syst.nbOfParticles();
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      output.obsInternalTemperature().currentValue() += 1/internalTemperature(syst(iOfParticle).internalEnergy());
    output.obsInternalTemperature().currentValue() /= syst.nbOfParticles();
  }
  
  void Dynamics::getThermo(Output& output) const
  {
    output.temperature() = 2 * output.kineticEnergy() / (output.dimension() * output.nbOfParticles());
    output.obsTotalEnergy().currentValue() = output.kineticEnergy() + output.potentialEnergy();
  }
  
  ///
  ///Computes the pressure from the kineticEnergy and the totalVirial, these fields must be updated
  void Dynamics::getPressure(Output& output) const
  {
    output.pressure() = (2 * output.kineticEnergy() + output.totalVirial()) / (output.dimension() * output.nbOfParticles() * pow(output.latticeParameter(), output.dimension()));
  }

}
