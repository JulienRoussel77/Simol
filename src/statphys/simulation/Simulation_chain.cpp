#include "Simulation.hpp"

namespace simol {

  void sampleSystem(BoundaryLangevin& dyna, BiChain& syst)
  {
    cout << "Initialization of the system...";cout.flush();
    double alpha, localTemp, localDist;
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
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

    for (int iOfIteration  =0; iOfIteration < dyna.nbOfThermalIterations(); ++iOfIteration)
      {
	syst.thermalize(dyna);
      }
    
    cout << "Done ! / Burning...";cout.flush();
    
    for (int iOfIteration  =0; iOfIteration < dyna.nbOfBurningIterations(); ++iOfIteration)
    {
      simulate(dyna, syst);
    }
    cout << "Done !" << endl;
  }
  
  template <>
  void computeOutput(BoundaryLangevin const& dyna, Chain const& syst, Output& output, int iOfIteration)
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
  
  void writeFinalOutput(BoundaryLangevin const& dyna, BiChain const& syst, Output& output)
  {
    output.finalDisplay(syst.configuration(), syst.externalForce());
    if (output.doComputeCorrelations())
      output.finalDisplayAutocorrelations();
    output.displayFinalFlow(dyna.temperature(), dyna.deltaTemperature());
  }


  void sampleSystem(BoundaryLangevin& dyna, TriChain& syst)
  {
    cout << "Initialization of the system...";cout.flush();
    double alpha, localTemp, localBending;
    //ofstream outTest("test.txt");
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
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

    for (int iOfIteration  =0; iOfIteration < dyna.nbOfThermalIterations(); ++iOfIteration)
    {
      syst.thermalize(dyna);
    }

    cout << "Done ! / Burning...";cout.flush();

    for (int iOfIteration  =0; iOfIteration < dyna.nbOfBurningIterations(); ++iOfIteration)
    {
      simulate(dyna, syst);
    }
    cout << "Done !" << endl;
  }
 
  ///
  ///Computes the quantities needed by the control variates (coefficients a, b, D) and {L \Phi}
  void updateAllControlVariates(const BoundaryLangevin& dyna, System const& syst, Output& output, int iOfIteration)
  {
    Vector<double> generatorOnBasis;
    
    //generatorOnBasis = generatorOn(dyna, syst, output.midFlowCV());
    generatorOnBasis = output.midFlowCV().generatorBoundarylangevin(syst.configuration(), dyna.betaLeft(), dyna.betaRight(), dyna.gamma());
    output.midFlowCV().update(output.energyMidFlow(), generatorOnBasis, syst.configuration(), iOfIteration);
    //generatorOnBasis = generatorOn(dyna, syst, output.sumFlowCV());
    generatorOnBasis = output.sumFlowCV().generatorBoundarylangevin(syst.configuration(), dyna.betaLeft(), dyna.betaRight(), dyna.gamma());
    output.sumFlowCV().update(output.energySumFlow(), generatorOnBasis, syst.configuration(), iOfIteration);
  }
  

  
  void writeOutput(BoundaryLangevin const& dyna, Chain const& syst, Output& output, int iOfIteration)
  {
    if (output.doOutput(iOfIteration))// && iOfIteration >= 100)
    {
      output.displayObservables(iOfIteration);
      output.displayChainPositions(syst.configuration(), iOfIteration);
      output.displayChainMomenta(syst.configuration(), iOfIteration);
      output.displayParticles(syst.configuration(), iOfIteration);
      
      output.midFlowCV_->display(output.outMidFlowCV(), iOfIteration * dyna.timeStep() );
      output.sumFlowCV_->display(output.outSumFlowCV(), iOfIteration * dyna.timeStep() );
    }
    
    if (output.doProfileOutput(iOfIteration))
    {
      output.displayProfile(iOfIteration);
      output.displayParticles(syst.configuration(), iOfIteration);
    }
  }
  
  void writeFinalOutput(BoundaryLangevin const& dyna, TriChain const& syst, Output& output)
  {
    output.finalDisplay(syst.configuration(), syst.externalForce());
    if (output.doComputeCorrelations())
      output.finalDisplayAutocorrelations();
    output.displayFinalFlow(dyna.temperature(), dyna.deltaTemperature(), dyna.tauBending(), dyna.xi());
  }
  
  void simulate(BoundaryLangevin& dyna, Chain& syst)
  {
    for (auto&& particle:syst.configuration())
      dyna.verletFirstPart(particle);

    syst.computeAllForces();

    for (auto&& particle:syst.configuration())
      dyna.verletSecondPart(particle);

    //for (int iOfParticle=0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
    //  dyna.updateAfter(getParticle(iOfParticle));

    dyna.updateOrsteinUhlenbeck(syst.getParticle(0), dyna.betaLeft());
    dyna.updateOrsteinUhlenbeck(syst.getParticle(syst.nbOfParticles() - 1), dyna.betaRight());

    if (dyna.doMomentaExchange())
      for (int iOfParticle=0; iOfParticle < syst.nbOfParticles()-1; iOfParticle++)
        dyna.updateMomentaExchange(syst.getParticle(iOfParticle), syst.getParticle(iOfParticle+1));
  }

}
