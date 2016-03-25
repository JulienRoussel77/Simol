#include "simulation.hpp"

using std::cout; 
using std::endl; 

namespace simol {
  
  void printClass(ParticleSystem& syst)
  {
    cout << "System = System" << endl;
  }
  
  void printClass(Isolated& syst)
  {
    cout << "System = Isolated" << endl;
  }
  
  void printClass(Langevin& syst)
  {
    cout << "Dynamics = Langevin" << endl;
  }
  
  void printClass(Dynamics& syst)
  {
    cout << "Dynamics = Dynamics" << endl;
  }
  
 Simulation::Simulation(Input& input):
  dimension_(input.dimension()),
  output_(input),
  //controlVariates_(nbOfReplicas_),
  rng_(input.seed(), input.dimension())
 {
   cout << "Simulation::Simulation(Input& input)" << endl;
   input.rng() = std::shared_ptr<RNG>(&rng_);
   output_.verbose() = 1;
   cout << "verbose = " << output_.verbose() << endl;
   system_ = createSystem(input); 
   cout << "Simulation : ";
   //system_->setRNG(&rng_);
   dynamics_ = createDynamics(input);
   output_.setControlVariates(input, system_->potential(), dynamics_->galerkin());
   
   cout << "----------" << endl;
   printClass(*system_);
   printClass(*dynamics_);
   system_->printName();
   dynamics_->printName();
   cout << "----------" << endl;
 }
 
 Simulation::~Simulation()
 {
   delete system_;
   delete dynamics_;
 }
 
 
 
 
  ParticleSystem* createSystem(Input  const& input)
  {
    cout << "createSystem(Input  const& input)" << endl;
    if (input.systemName() == "Isolated")
      return new Isolated(input);
    else if (input.systemName() == "BiChain")
      return new BiChain(input);
    else if (input.systemName() == "TriChain")
      return new TriChain(input);
    else
      std::cout << input.systemName() << " is not a valid system !" << std::endl;
    
    return 0;
  }

	///
	///Calls the proper constructor of Dynamics and returns a pointer
  Dynamics* createDynamics(Input const& input)
  {
		//cout << "createDynamics(Input const& input)" << endl;
    if (input.dynamicsName() == "Hamiltonian")
      return new Hamiltonian(input);
    else if (input.dynamicsName() == "Langevin")
      return new Langevin(input);
    else if (input.dynamicsName() == "BoundaryLangevin")
      return new BoundaryLangevin(input);
    else if (input.dynamicsName() == "Overdamped")
      return new Overdamped(input);
    else if (input.dynamicsName() == "DPDE")
      return new DPDE(input);
    else
      std::cout << input.dynamicsName() << " is not a valid dynamics !" << std::endl;
    return 0;
  }
  
  
 

 
  ///
  ///Initializes all the momenta, assuming an affine temperature profile
  void initializeMomenta(const Dynamics& dyna, ParticleSystem& syst)
  {}
  
  void simulate(Dynamics& dyna, ParticleSystem& syst)
  {
    for (auto&& particle : syst.configuration())
      dyna.updateBefore(particle);
    
    syst.computeAllForces(dyna);
    
    for (auto&& particle : syst.configuration())
      dyna.updateAfter(particle);
  }
  
  void initializeSystem(Dynamics& dyna, ParticleSystem& syst) {throw std::invalid_argument("initializeSystem : Function undefined");}
 
  ///
  ///Computes the quantities needed by the control variates (coefficients a, b, D) and {L \Phi}
  void updateAllControlVariates(const Dynamics& dyna, const ParticleSystem& syst, Output& output, size_t iOfIteration)
  {
    Vector<double> q = syst.getParticle(0).position();
    Vector<double> p = syst.getParticle(0).momentum();
    Vector<double> qEnd = syst.getParticle(syst.nbOfParticles()-1).position();
    Vector<double> generatorOnBasis;
    generatorOnBasis = generatorOn(dyna, syst, output.velocityCV());
    output.velocityCV().update(p(0), generatorOnBasis, syst.configuration(), iOfIteration);
    generatorOnBasis = generatorOn(dyna, syst, output.forceCV());
    output.forceCV().update(syst.force(q)(0), generatorOnBasis, syst.configuration(), iOfIteration);
    generatorOnBasis = generatorOn(dyna, syst, output.lengthCV());
    output.lengthCV().update(qEnd(0), generatorOnBasis, syst.configuration(), iOfIteration);
    generatorOnBasis = generatorOn(dyna, syst, output.midFlowCV());
    output.midFlowCV().update(output.energyMidFlow(), generatorOnBasis, syst.configuration(), iOfIteration);
    generatorOnBasis = generatorOn(dyna, syst, output.sumFlowCV());
    output.sumFlowCV().update(output.energySumFlow(), generatorOnBasis, syst.configuration(), iOfIteration);
  }
  
  ///
  ///Applies the generator of this dynamics to the basis functions of the CV
  Vector<double> generatorOn(const Dynamics& dyna, const ParticleSystem& syst, const ControlVariate& controlVariate)
  {
    throw std::invalid_argument("GeneratorOn not defined in the general case");  
  }
  
  void computeOutput(const Dynamics& dyna, const ParticleSystem& syst, Output& output, size_t iOfIteration)
  {
    if (output.verbose() > 0)
    {
      output.kineticEnergy() = 0;
      output.potentialEnergy() = 0;
      //Calcul de la température et de l'énergie
      for (const auto& particle : syst.configuration())
      {
        //Particle& particle = syst.getParticle(iOfParticle);
        output.kineticEnergy() += particle.kineticEnergy();
        output.potentialEnergy() += particle.potentialEnergy();
      }
    }
    // In the case of the trichain we add the potential of the wall interaction
    output.potentialEnergy() += syst.boundaryPotEnergy();
    syst.computeProfile(output, dyna, iOfIteration);
    updateAllControlVariates(dyna, syst, output, iOfIteration);
  }
  
  //##################### ISOLATED ######################
  
  void initializeSystem(Dynamics& dyna, Isolated& syst)
  {
    syst.getParticle(0).momentum() = syst.drawMomentum(dyna.beta(), syst.getParticle(0).mass());
    syst.getParticle(0).position(0) = syst.drawPotLaw(dyna.beta()); 
  }
  
    void initializeSystem(Langevin& dyna, Isolated& syst)
  {
    cout << "youpi" << endl;
  }
  
  //##################### HAMILTONIAN ####################"
  
    ///
  ///Applies the generator of this dynamics to the basis functions of the CV
  Vector<double> generatorOn(const Hamiltonian& dyna, const ParticleSystem& syst, ControlVariate const& controlVariate)
  {
    Vector<double> result = Vector<double>::Zero(controlVariate.nbOfFunctions());
    for (size_t iOfFunction=0; iOfFunction < controlVariate.nbOfFunctions(); iOfFunction++)
      for (size_t iOfParticle=0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
        result(iOfFunction) += syst.getParticle(iOfParticle).momentum().dot(controlVariate.gradientQ(syst.configuration(), iOfParticle, iOfFunction))
        + syst.force(syst.getParticle(iOfParticle).position()).dot(controlVariate.gradientP(syst.configuration(), iOfParticle, iOfFunction))
        + dyna.externalForce().dot(controlVariate.gradientP(syst.configuration(), iOfParticle, iOfFunction));  
    return result;      
  }
  
  //################# LANGEVIN ###########################
  
  ///Applies the generator of this dynamics to the basis functions of the CV
  ///Evaluate at the current state of "conifguration"
  Vector<double> generatorOn(const Langevin& dyna, const ParticleSystem& syst, const ControlVariate& controlVariate)
  {
    Vector<double> result = Vector<double>::Zero(controlVariate.nbOfFunctions());
    for (size_t iOfFunction=0; iOfFunction < controlVariate.nbOfFunctions(); iOfFunction++)
      for (size_t iOfParticle=0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      {
        result(iOfFunction) += dot(syst.getParticle(iOfParticle).momentum(), controlVariate.gradientQ(syst.configuration(), iOfParticle, iOfFunction))
          + dot(syst.getParticle(iOfParticle).force(), controlVariate.gradientP(syst.configuration(), iOfParticle, iOfFunction))
          + dyna.gamma() * (- dot(syst.getParticle(iOfParticle).momentum(), controlVariate.gradientP(syst.configuration(), iOfParticle, iOfFunction))
              + controlVariate.laplacianP(syst.configuration(), iOfParticle, iOfFunction) / dyna.beta() );
      }
    return result;
  }
  
    ///
  ///Computes the quantities needed by the control variates (coefficients a, b, D) and {L \Phi}
  void updateAllControlVariates(const Langevin& dyna, const ParticleSystem& syst, Output& output, size_t iOfIteration)
  {
    //cout << "Langevin::updateAllControlVariates(Output& output, vector<Particle> const& configuration, size_t iOfIteration)" << endl;
    Vector<double> q = syst.getParticle(0).position();
    Vector<double> p = syst.getParticle(0).momentum();
    Vector<double> generatorOnBasis;
    generatorOnBasis = generatorOn(dyna, syst, output.velocityCV());
    output.velocityCV().update(p(0), generatorOnBasis, syst.configuration(), iOfIteration);
    if (output.doOutput(iOfIteration))
      output.displayGeneratorOnBasis(output.outVelocitiesGenerator_, syst.configuration(), output.velocityCV(), iOfIteration*dyna.timeStep());
    
    /*generatorOnBasis = generatorOn(output.forceCV(), configuration);
    output.forceCV().update(potential_->derivative(q)(0), generatorOnBasis, configuration, iOfIteration);
    generatorOnBasis = generatorOn(output.lengthCV(), configuration);
    output.lengthCV().update(q(0), generatorOnBasis, configuration, iOfIteration);*/
    //cout << "end updateAllControlVariates" << endl;
  }
  
  //################### OVERDAMPED ########################
  
  ///Applies the generator of this dynamics to the basis functions of the CV
  ///Evaluate at the current state of "conifguration"
  Vector<double> generatorOn(const Overdamped& dyna, const ParticleSystem& syst, const ControlVariate& controlVariate)
  {
    Vector<double> result = Vector<double>::Zero(controlVariate.nbOfFunctions());
    for (size_t iOfFunction=0; iOfFunction < controlVariate.nbOfFunctions(); iOfFunction++)
      for (size_t iOfParticle=0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
        result(iOfFunction) += controlVariate.laplacianQ(syst.configuration(), iOfParticle, iOfFunction) / dyna.beta()
          + dot(syst.force(syst.getParticle(iOfParticle).position()), controlVariate.gradientQ(syst.configuration(), iOfParticle, iOfFunction));
    return result;
  }
  
  //################### CHAINS #############################
  
  ///
  ///Initializes all the momenta, assuming an affine temperature profile
  void initializeMomenta(const BoundaryLangevin& dyna, Chain& syst)
  {
    for (size_t iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
    {
      Particle& particle = syst.getParticle(iOfParticle);
      double tempi_ = dyna.temperatureLeft() + (iOfParticle + .5) * (dyna.temperatureRight() - dyna.temperatureLeft()) / syst.nbOfParticles();
      particle.momentum() = sqrt(tempi_ / particle.mass()) * syst.rng()->gaussian();
    }
  }
  
  void initializeSystem(BoundaryLangevin& dyna, BiChain& syst)
  {
    cout << "Initialization of the system...";cout.flush();
    double alpha, localTemp, localDist;
    //ofstream outTest("test.txt");
    for (size_t iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
    {
      alpha = iOfParticle / (double) syst.nbOfParticles();     
      localTemp = (1-alpha) * dyna.temperatureLeft() + alpha * dyna.temperatureRight();
      //localBending = dyna.computeMeanPotLaw(1/localTemp);
      localDist = syst.drawPotLaw(1/localTemp);
      //outTest << iOfParticle << " " << localBending << endl;
      //cout << "bending = " << localBending << " / mean = " << dyna.computeMeanPotLaw(1/localTemp) << endl;
      
      double prevPosition = (iOfParticle>0)?syst.getParticle(iOfParticle-1).position(0):0;
      
      
      syst.getParticle(iOfParticle).position(0) = prevPosition + localDist;
      //cout << -position2 << " + 2 * " << position1 << " + " << localBending << " = " << syst.getParticle(iOfParticle).position(0) << endl;
      syst.getParticle(iOfParticle).momentum() = syst.drawMomentum(1/localTemp, syst.getParticle(iOfParticle).mass());
    
      dyna.initializeCountdown(syst.getParticle(iOfParticle));
    }
    
    cout << "Done ! / Thermalization...";cout.flush();
    
    for (size_t iOfIteration  =0; iOfIteration < dyna.nbOfThermalIterations(); ++iOfIteration)
    {
      syst.thermalize(dyna);
    }
    
    cout << "Done ! / Burning...";cout.flush();
    
    for (size_t iOfIteration  =0; iOfIteration < dyna.nbOfBurningIterations(); ++iOfIteration)
    {
      simulate(dyna, syst);
    }
    cout << "Done !" << endl;
  }
  
  void initializeSystem(BoundaryLangevin& dyna, TriChain& syst)
  {
    cout << "Initialization of the system...";cout.flush();
    double alpha, localTemp, localBending;
    //ofstream outTest("test.txt");
    for (size_t iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
    {
      alpha = iOfParticle / (double) syst.nbOfParticles();     
      localTemp = (1-alpha) * dyna.temperatureLeft() + alpha * dyna.temperatureRight();
      //localBending = dyna.computeMeanPotLaw(1/localTemp);
      localBending = syst.drawPotLaw(1/localTemp);
      //outTest << iOfParticle << " " << localBending << endl;
      //cout << "bending = " << localBending << " / mean = " << dyna.computeMeanPotLaw(1/localTemp) << endl;
      
      double position1 = (iOfParticle>0)?syst.getParticle(iOfParticle-1).position(0):0;
      double position2 = (iOfParticle>1)?syst.getParticle(iOfParticle-2).position(0):0;
      
      
      syst.getParticle(iOfParticle).position(0) = -position2 + 2 * position1 + localBending;
      //cout << -position2 << " + 2 * " << position1 << " + " << localBending << " = " << syst.getParticle(iOfParticle).position(0) << endl;
      syst.getParticle(iOfParticle).momentum() = syst.drawMomentum(1/localTemp, syst.getParticle(iOfParticle).mass());
      
      dyna.initializeCountdown(syst.getParticle(iOfParticle));
    }
    
    cout << "Done ! / Thermalization...";cout.flush();
    
    for (size_t iOfIteration  =0; iOfIteration < dyna.nbOfThermalIterations(); ++iOfIteration)
    {
      syst.thermalize(dyna);
    }
    
    cout << "Done ! / Burning...";cout.flush();
    
    for (size_t iOfIteration  =0; iOfIteration < dyna.nbOfBurningIterations(); ++iOfIteration)
    {
      simulate(dyna, syst);
    }
    cout << "Done !" << endl;
  }
  
  
  void simulate(BoundaryLangevin& dyna, Chain& syst)
  {
    //for (auto&& particle : configuration_)
      //dyna.updateBefore(particle);
    //for (size_t iOfParticle=0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
    //  dyna.updateBefore(getParticle(iOfParticle));
    
    for (auto&& particle:syst.configuration())
      dyna.updateBefore(particle);
    
    syst.computeAllForces(dyna);
    
    for (auto&& particle:syst.configuration())
      dyna.updateAfter(particle);
    
    //for (size_t iOfParticle=0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
    //  dyna.updateAfter(getParticle(iOfParticle));
    
    dyna.updateOrsteinUhlenbeck(syst.getParticle(0), dyna.betaLeft());
    dyna.updateOrsteinUhlenbeck(syst.getParticle(syst.nbOfParticles() - 1), dyna.betaRight());
    
    if (dyna.doMomentaExchange())
      for (size_t iOfParticle=0; iOfParticle < syst.nbOfParticles()-1; iOfParticle++)
        dyna.updateMomentaExchange(syst.getParticle(iOfParticle), syst.getParticle(iOfParticle+1));
  }
  
  
    ///Applies the generator of this dynamics to the basis functions of the CV
  ///Evaluate at the current state of "conifguration"
  Vector<double> generatorOn(const BoundaryLangevin& dyna, const ParticleSystem& syst, const ControlVariate& controlVariate)
  {
    size_t nbOfParticles = syst.nbOfParticles();
    Vector<double> result = Vector<double>::Zero(controlVariate.nbOfFunctions());
    for (size_t iOfFunction=0; iOfFunction < controlVariate.nbOfFunctions(); iOfFunction++)
    {
      for (size_t iOfParticle=0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      result(iOfFunction) += dot(syst.getParticle(iOfParticle).momentum(), controlVariate.gradientQ(syst.configuration(), iOfParticle, iOfFunction))
        + dot(syst.getParticle(iOfParticle).force(), controlVariate.gradientP(syst.configuration(), iOfParticle, iOfFunction));
        //if(false)       
      result(iOfFunction) += dyna.gamma() * (- dot(syst.getParticle(0).momentum(), controlVariate.gradientP(syst.configuration(), 0, iOfFunction))
      + controlVariate.laplacianP(syst.configuration(), 0, iOfFunction) / dyna.betaLeft()
      - dot(syst.getParticle(nbOfParticles-1).momentum(), controlVariate.gradientP(syst.configuration(), nbOfParticles-1, iOfFunction))
      + controlVariate.laplacianP(syst.configuration(), nbOfParticles-1, iOfFunction) / dyna.betaRight());
    }
    return result;   
  }

  //-------------- algorithmes pour DPDE -------------
  
  void simulate(DPDE& dyna, ParticleSystem& syst)
  {
    //-- Verlet part -- TO DO : a completer...
    //for (auto&& particle : syst.configuration())
    //dyna.verletFirstPart(particle);
    syst.computeAllForces(dyna);
    //for (auto&& particle : syst.configuration())
    //dyna.verletSecondPart(particle);
    //-- fluctuation/dissipation --
    //for (auto&& particle : syst.configuration())
    for (int i=0; i < syst.nbOfParticles(); i++)  // version explicite de la ligne cachee ci dessus, utile pour faire des boucles sur les couples !
      dyna.energyReinjection(syst.getParticle(i));  // integration de p avec gamma fixe + reinjection
  }
  
  void initializeSystem(DPDE& dyna, Isolated& syst)
  {
    syst.getParticle(0).momentum() = syst.drawMomentum(dyna.beta(), syst.getParticle(0).mass());
    syst.getParticle(0).position(0) = 0;
    syst.getParticle(0).internalEnergy() = 1;  // TO DO : il faudra ici tirer selon la bonne loi ?
  }

  void computeOutput(const DPDE& dyna, const Isolated& syst, Output& output, size_t iOfIteration)
  {
    output.kineticEnergy() = syst.getParticle(0).kineticEnergy();
    output.potentialEnergy() = syst.getParticle(0).potentialEnergy();
    output.internalEnergy() = syst.getParticle(0).internalEnergy();
  }

  
    //------------------------------ MAIN FUNCTION ------------------------------------
  template<typename D, typename S>
  void launchSimu(D& dyna, S& syst, Output& output)
  {
    //dyna.printName();
    //syst.printName();
    printClass(dyna);
    printClass(syst);
    cout << "Estimated time : " << 3.5 * syst.nbOfParticles()/1024. * dyna.nbOfIterations() / 1e6 << " hours" << endl;
        
    //---- initialization (including burn-in) -----
    initializeSystem(dyna, syst);
    syst.computeAllForces(dyna);  // TO DO : a mettre dans la fonction d'initialisation...
    
    //---- actual iterations -----
    for (size_t iOfIteration  =0; iOfIteration < dyna.nbOfIterations(); ++iOfIteration)
      {
        //--- display progress every time 10% of simulation elapsed ---
        if ((10*iOfIteration) % dyna.nbOfIterations() == 0) 
          cout << "---- Run " << (100 * iOfIteration) / dyna.nbOfIterations() << " % completed ----" << endl;
      
        //--- write outputs if required ----
        computeOutput(dyna, syst, output, iOfIteration);
        syst.writeOutput(output, iOfIteration);

        //---- update the system by the numerical integration ---
        simulate(dyna, syst);
      }
    
    //--- write final outputs ----
    syst.computeFinalOutput(output, dyna);
    syst.writeFinalOutput(output, dyna);
  }
  

  //------------ fonction de lancement ----------------

  void Simulation::launch()  
  {
    cout << "Simulation::launch" << endl;
    //dynamics_->printName();
    //system_->printName();
    printClass(*dynamics_);
    printClass(*system_);
    launchSimu(*dynamics_, *system_, output_);
  }
  


  
  
  
  
}
