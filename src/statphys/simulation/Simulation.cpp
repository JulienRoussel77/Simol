#include "Simulation.hpp"

namespace simol {
 
  
  //void writeFinalOutput(Dynamics const& dyna, System const& syst, Output& output);
  
  //###### Global ######
  

  void sampleSystem(Dynamics& dyna, System& syst)
  {
    dyna.printName();
    syst.printName();
    throw std::invalid_argument("sampleSystem : Function undefined");
  }
  
  ///
  ///Computes the quantities needed by the control variates (coefficients a, b, D) and {L \Phi}
  void updateAllControlVariates(Dynamics const& /*dyna*/, System const& /*syst*/, Output& /*output*/, size_t /*iOfIteration*/)
  {throw std::invalid_argument("updateAllControlVariates: Function undefined");}

    /*Vector<double> q = syst.getParticle(0).position();
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
  }*/
    
      ///
  ///Applies the generator of this dynamics to the basis functions of the CV
  Vector<double> generatorOn(Dynamics const& /*dyna*/, System const& /*syst*/, const ControlVariate& /*controlVariate*/)
  {
    throw std::invalid_argument("GeneratorOn not defined in the general case");
  }

  void computeOutput(Dynamics const& dyna, System const& syst, Output& output, size_t iOfIteration)
  {
    output.kineticEnergy() = 0;
    output.potentialEnergy() = 0;
    //Calcul de la température et de l'énergie
    for (const auto& particle : syst.configuration())
    {
      output.kineticEnergy() += particle.kineticEnergy();
      output.potentialEnergy() += particle.potentialEnergy();
    }
    // In the case of the trichain we add the potential of the wall interaction
    output.potentialEnergy() += syst.boundaryPotEnergy();
    syst.computeProfile(output, dyna, iOfIteration);
    updateAllControlVariates(dyna, syst, output, iOfIteration);
  }

  // TO DO ou pas: cette fonction ne sert pas a grand chose... ce test peut (DOIT) etre fait dans output.cpp !!
  void writeOutput(Dynamics const& /*dyna*/, System const& /*syst*/, Output& /*output*/, size_t /*iOfIteration*/)
  {
    throw std::invalid_argument("writeOutput not defined in the general case");
  }

  void writeFinalOutput(Dynamics const& dyna, System const& syst, Output& output)
  {
    output.finalDisplay(syst.configuration(), dyna.externalForce());
    if (output.doComputeCorrelations())
      output.finalDisplayAutocorrelations();
  }

  void simulate(Dynamics& dyna, System& syst)
  {
    for (auto&& particle : syst.configuration())
      dyna.updateBefore(particle);

    syst.computeAllForces(dyna);

    for (auto&& particle : syst.configuration())
      dyna.updateAfter(particle);
  }

  
  //###### ISOLATED #####
  
  void sampleSystem(Hamiltonian& /*dyna*/, Isolated& /*syst*/)
  {}
  
  void sampleSystem(UniformStochasticDynamics& dyna, Isolated& syst)
  {
    syst.getParticle(0).momentum() = syst.drawMomentum(dyna.beta(), syst.getParticle(0).mass());
    syst.getParticle(0).position(0) = syst.drawPotLaw(dyna.beta());
  }
  

  void sampleSystem(Overdamped& dyna, Isolated& syst)
  {
    syst.getParticle(0).position(0) = syst.drawPotLaw(dyna.beta());
  }
  
  /*void writeFinalOutput(Dynamics const& dyna, Isolated const& syst, Output& output)
  {
    output.finalDisplay(syst.configuration(), dyna.externalForce());
    if (output.doComputeCorrelations())
      output.finalDisplayAutocorrelations();
    output.displayFinalVelocity(dyna.temperature(), dyna.externalForce(0), output.velocityCV_->nbOfFourier(), output.velocityCV_->nbOfHermite());
  }*/
  
  void computeOutput(Dynamics const& dyna, const Isolated& syst, Output& output, size_t iOfIteration)
  {
    //cout << "computeOutput(D const& dyna, const Isolated& syst, Output& output, size_t iOfIteration)" << endl;
    //Calcul de la température et de l'énergie
    output.kineticEnergy() = syst.getParticle(0).kineticEnergy();
    output.potentialEnergy() = syst.getParticle(0).potentialEnergy();
    updateAllControlVariates(dyna, syst, output, iOfIteration);
  }
  
  void writeFinalOutput(Hamiltonian const& dyna, Isolated const& syst, Output& output)
  {
    output.finalDisplay(syst.configuration(), dyna.externalForce());
    if (output.doComputeCorrelations())
      output.finalDisplayAutocorrelations();
    output.displayFinalVelocity(0, dyna.externalForce(0), output.velocityCV_->nbOfFourier(), output.velocityCV_->nbOfHermite());
  }
  
  
  //###### CHAIN ######
  
  void sampleSystem(BoundaryLangevin& dyna, BiChain& syst)
  {
    cout << "Initialization of the system...";cout.flush();
    double alpha, localTemp, localDist;
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
  
  void writeFinalOutput(BoundaryLangevin const& dyna, BiChain const& syst, Output& output)
  {
    output.finalDisplay(syst.configuration(), dyna.externalForce());
    if (output.doComputeCorrelations())
      output.finalDisplayAutocorrelations();
    output.displayFinalFlow(dyna.temperature(), dyna.deltaTemperature());
  }


  void sampleSystem(BoundaryLangevin& dyna, TriChain& syst)
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
 
  ///Applies the generator of this dynamics to the basis functions of the CV
  ///Evaluate at the current state of "conifguration" 
  Vector<double> generatorOn(const BoundaryLangevin& dyna, System const& syst, const ControlVariate& controlVariate)
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
  
  ///
  ///Computes the quantities needed by the control variates (coefficients a, b, D) and {L \Phi}
  void updateAllControlVariates(const BoundaryLangevin& dyna, System const& syst, Output& output, size_t iOfIteration)
  {
    Vector<double> generatorOnBasis;
    
    generatorOnBasis = generatorOn(dyna, syst, output.midFlowCV());
    output.midFlowCV().update(output.energyMidFlow(), generatorOnBasis, syst.configuration(), iOfIteration);
    generatorOnBasis = generatorOn(dyna, syst, output.sumFlowCV());
    output.sumFlowCV().update(output.energySumFlow(), generatorOnBasis, syst.configuration(), iOfIteration);
  }
  
  void writeOutput(BoundaryLangevin const& dyna, System const& syst, Output& output, size_t iOfIteration)
  {
    if (output.doOutput(iOfIteration))// && iOfIteration >= 100)
    {
      output.displayObservables(iOfIteration);
      output.displayChainPositions(syst.configuration(), iOfIteration);
      output.displayChainMomenta(syst.configuration(), iOfIteration);
      output.displayParticles(syst.configuration(), iOfIteration);
      
      output.midFlowCV_->display(output.outMidFlowCV_, iOfIteration * dyna.timeStep() );
      output.sumFlowCV_->display(output.outSumFlowCV_, iOfIteration * dyna.timeStep() );
    }
    
    if (output.doProfileOutput(iOfIteration))
    {
      output.displayProfile(iOfIteration);
      output.displayParticles(syst.configuration(), iOfIteration);
    }
  }
  
  void writeFinalOutput(BoundaryLangevin const& dyna, TriChain const& syst, Output& output)
  {
    output.finalDisplay(syst.configuration(), dyna.externalForce());
    if (output.doComputeCorrelations())
      output.finalDisplayAutocorrelations();
    output.displayFinalFlow(dyna.temperature(), dyna.deltaTemperature(), dyna.tauBending(), dyna.xi());
  }
  
  void simulate(BoundaryLangevin& dyna, Chain& syst)
  {
    for (auto&& particle:syst.configuration())
      dyna.verletFirstPart(particle);

    syst.computeAllForces(dyna);

    for (auto&& particle:syst.configuration())
      dyna.verletSecondPart(particle);

    //for (size_t iOfParticle=0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
    //  dyna.updateAfter(getParticle(iOfParticle));

    dyna.updateOrsteinUhlenbeck(syst.getParticle(0), dyna.betaLeft());
    dyna.updateOrsteinUhlenbeck(syst.getParticle(syst.nbOfParticles() - 1), dyna.betaRight());

    if (dyna.doMomentaExchange())
      for (size_t iOfParticle=0; iOfParticle < syst.nbOfParticles()-1; iOfParticle++)
        dyna.updateMomentaExchange(syst.getParticle(iOfParticle), syst.getParticle(iOfParticle+1));
  }
  
  //------------------- NBody ---------------------
  
    void sampleSystem(Dynamics& dyna, NBody& syst)
  {
    int Dim = syst.dimension();   // PAS SUPER, MAIS SINON PBM DE TYPE POUR COMPARAISON ?
    int NbPartDim = syst.nbOfParticlesPerDimension(); 
    double latticeSize = syst.latticeParameter();
    //-- initialization of the momenta according to a Gaussian distribution --
    for (size_t i = 0; i < syst.nbOfParticles(); i++)
      syst.getParticle(i).momentum() = syst.drawMomentum(1, syst.getParticle(i).mass());   // TO DO : introduce beta in Hamiltonian...
    //-- initialization on a cubic lattice --
    if (Dim == 2) 
      {
	for (int i = 0; i < NbPartDim; i++)
	  for (int j = 0; j < NbPartDim; j++)
	    {
	      syst.getParticle(i*NbPartDim+j).position(0) = i*latticeSize;
	      syst.getParticle(i*NbPartDim+j).position(1) = j*latticeSize;
	    }
      }
    else if (Dim == 3) 
      {
	int NbPartDim2 = NbPartDim*NbPartDim;
	for (int i = 0; i < NbPartDim; i++)
	  for (int j = 0; j < NbPartDim; j++)
	    for (int k = 0; k < NbPartDim; k++)
	    {
	      syst.getParticle(i*NbPartDim2+j*NbPartDim+k).position(0) = i*latticeSize;
	      syst.getParticle(i*NbPartDim2+j*NbPartDim+k).position(1) = j*latticeSize;
	      syst.getParticle(i*NbPartDim2+j*NbPartDim+k).position(2) = k*latticeSize;
	    }
      }
    else 
      {
	throw std::invalid_argument("sampleSystem: Bad dimension, should be 2 or 3");
      }
    //cout << "    VERIFICATION : " << syst.nbOfParticlesPerDimension() << endl;
  }
  
  void simulate(Hamiltonian& dyna, NBody& syst)
  {
    for (auto&& particle : syst.configuration())
      dyna.verletFirstPart(particle);
    syst.computeAllForces(dyna);
    for (auto&& particle : syst.configuration())
      dyna.verletSecondPart(particle);
  }
  
  void computeOutput(Hamiltonian const& /*dyna*/, NBody const& syst, Output& output, size_t /*iOfIteration*/)
  {
    output.kineticEnergy() = 0;
    output.potentialEnergy() = 0;
    output.totalVirial() = 0;
    //Calcul de la température et de l'énergie
    for (const auto& particle : syst.configuration())
      {
      output.kineticEnergy() += particle.kineticEnergy();
      output.potentialEnergy() += particle.potentialEnergy();	      
      output.totalVirial() += particle.virial();
      }
  }
  
  //--- CONFLIT DE TEMPLETAGE ICI AUSSI : entre dynamics et system... on ne peut pas preciser que le systeme ?! ---
  void writeOutput(Hamiltonian const& /*dyna*/, NBody const& syst, Output& output, size_t iOfIteration)
  {
    if (output.doOutput(iOfIteration))
      output.displayObservables(iOfIteration);
    if (output.doProfileOutput(iOfIteration))
      output.displayParticlesXMakeMol(syst.configuration(), iOfIteration, syst.latticeParameter()*syst.nbOfParticlesPerDimension());        
  }
  
  void writeFinalOutput(Hamiltonian const& dyna, NBody const& syst, Output& output)
  {
    //output.finalDisplay(syst.configuration(), dyna.externalForce());
  }
  
  // ------------------------ Classified by Dynamics ---------------------------------
  
    //##################### HAMILTONIAN ####################"

  ///
  ///Applies the generator of this dynamics to the basis functions of the CV
  Vector<double> generatorOn(const Hamiltonian& dyna, System const& syst, ControlVariate const& controlVariate)
  {
    Vector<double> result = Vector<double>::Zero(controlVariate.nbOfFunctions());
    for (size_t iOfFunction=0; iOfFunction < controlVariate.nbOfFunctions(); iOfFunction++)
      for (size_t iOfParticle=0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
        result(iOfFunction) += syst.getParticle(iOfParticle).momentum().dot(controlVariate.gradientQ(syst.configuration(), iOfParticle, iOfFunction))
        + syst.force(syst.getParticle(iOfParticle).position()).dot(controlVariate.gradientP(syst.configuration(), iOfParticle, iOfFunction))
        + dyna.externalForce().dot(controlVariate.gradientP(syst.configuration(), iOfParticle, iOfFunction));
    return result;
  }

  void updateAllControlVariates(const Hamiltonian& /*dyna*/, System const& /*syst*/, Output& /*output*/, size_t /*iOfIteration*/)
  {}

  
  void writeOutput(Hamiltonian const& /*dyna*/, System const& syst, Output& output, size_t iOfIteration)
  {
    if (output.doOutput(iOfIteration))// && iOfIteration >= 100)
      output.displayObservables(iOfIteration);
    
    if (output.doProfileOutput(iOfIteration))
      output.displayParticles(syst.configuration(), iOfIteration);
  }
  
  // ###### OVERDAMPED ######
  
    ///Applies the generator of this dynamics to the basis functions of the CV
  ///Evaluate at the current state of "conifguration"
  Vector<double> generatorOn(const Overdamped& dyna, System const& syst, const ControlVariate& controlVariate)
  {
    Vector<double> result = Vector<double>::Zero(controlVariate.nbOfFunctions());
    for (size_t iOfFunction=0; iOfFunction < controlVariate.nbOfFunctions(); iOfFunction++)
      for (size_t iOfParticle=0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
        result(iOfFunction) += controlVariate.laplacianQ(syst.configuration(), iOfParticle, iOfFunction) / dyna.beta()
          + dot(syst.force(syst.getParticle(iOfParticle).position()), controlVariate.gradientQ(syst.configuration(), iOfParticle, iOfFunction));
    return result;
  }
  
  // ###### LANGEVIN ######
  
    ///Applies the generator of this dynamics to the basis functions of the CV
  ///Evaluate at the current state of "conifguration"
  Vector<double> generatorOn(const Langevin& dyna, System const& syst, const ControlVariate& controlVariate)
  {
    //cout << "generatorOn(const Langevin& dyna, S const& syst, const ControlVariate& controlVariate)" << endl;
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
  void updateAllControlVariates(const Langevin& dyna, System const& syst, Output& output, size_t iOfIteration)
  {
    //cout << "updateAllControlVariates(const Langevin& dyna, S const& syst, Output& output, size_t iOfIteration)" << endl;
    Vector<double> q = syst.getParticle(0).position();
    Vector<double> p = syst.getParticle(0).momentum();
    Vector<double> generatorOnBasis;
    generatorOnBasis = generatorOn(dyna, syst, output.velocityCV());
    output.velocityCV().update(p(0), generatorOnBasis, syst.configuration(), iOfIteration);
    if (output.doOutput(iOfIteration))
      output.displayGeneratorOnBasis(output.outVelocitiesGenerator_, syst.configuration(), output.velocityCV(), iOfIteration*dyna.timeStep());
  }
  
  void writeOutput(Langevin const& /*dyna*/, System const& syst, Output& output, size_t iOfIteration)
  {
    if (output.doOutput(iOfIteration))// && iOfIteration >= 100)
      output.displayObservables(iOfIteration);
    
    if (output.doProfileOutput(iOfIteration))
      output.displayParticles(syst.configuration(), iOfIteration);
  }
  
  // ###### DPDE ######
  
  void sampleSystem(DPDE& dyna, Isolated& syst)
  {
    syst.getParticle(0).momentum() = syst.drawMomentum(dyna.beta(), syst.getParticle(0).mass());
    syst.getParticle(0).position(0) = 0;
    syst.getParticle(0).internalEnergy() = 1;  // TO DO : il faudra ici tirer selon la bonne loi ?
  }
  
  void computeOutput(const DPDE& /*dyna*/, const Isolated& syst, Output& output, size_t /*iOfIteration*/)
  {
    output.kineticEnergy() = syst.getParticle(0).kineticEnergy();
    output.potentialEnergy() = syst.getParticle(0).potentialEnergy();
    output.internalEnergy() = syst.getParticle(0).internalEnergy();
  }
  
  void writeOutput(DPDE const& /*dyna*/, System const& syst, Output& output, size_t iOfIteration)
  {
    if (output.doOutput(iOfIteration))
      output.displayObservablesDPDE(syst.configuration(), iOfIteration);
  }
  
  void simulate(DPDE& dyna, System& syst)
  {
    for (auto&& particle : syst.configuration())
      dyna.verletFirstPart(particle);
    syst.computeAllForces(dyna);
    for (auto&& particle : syst.configuration())
      dyna.verletSecondPart(particle);
    //-- fluctuation/dissipation --
    for (size_t i=0; i < syst.nbOfParticles(); i++)  // version explicite de la ligne cachee ci dessus, utile pour faire des boucles sur les couples !
      dyna.energyReinjection(syst.getParticle(i));  // integration de p avec gamma fixe + reinjection
  }
  
  

  
  
  

  
  
  
}
