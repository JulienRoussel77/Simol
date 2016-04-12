#include "Simulation.hpp"

namespace simol {
  
  void sampleMomenta(BoundaryLangevin& dyna, Chain& syst)
  {
    cout << "Sampling the momenta...";cout.flush();
    double alpha, localTemp;
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      {
        alpha = iOfParticle / (double) syst.nbOfParticles();
        localTemp = (1-alpha) * dyna.temperatureLeft() + alpha * dyna.temperatureRight();
        syst.getParticle(iOfParticle).momentum() = syst.drawMomentum(1/localTemp, syst.getParticle(iOfParticle).mass());

        dyna.initializeCountdown(syst.getParticle(iOfParticle));
      }
      cout << "Done." << endl;
  }
  
  void samplePositions(BoundaryLangevin& dyna, Chain& syst)
  {
    throw std::invalid_argument("samplePositions(BoundaryLangevin& dyna, Chain& syst) not defined");
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
    
    if (output.doLongOutput(iOfIteration))
    {
      output.displayProfile(iOfIteration);
      output.displayParticles(syst.configuration(), iOfIteration);
    }
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
  
  
  //-------------- BiChain -----------------
  
  void samplePositions(BoundaryLangevin& dyna, BiChain& syst)
  {
    cout << "Sampling the positions..." ;cout.flush();
    double alpha, localTemp, localDist;
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      {
        alpha = iOfParticle / (double) syst.nbOfParticles();
        localTemp = (1-alpha) * dyna.temperatureLeft() + alpha * dyna.temperatureRight();
        localDist = syst.drawPotLaw(1/localTemp);
        double prevPosition = (iOfParticle>0)?syst.getParticle(iOfParticle-1).position(0):0;
      
        syst.getParticle(iOfParticle).position(0) = prevPosition + localDist;
      }
  }
  
  void writeFinalOutput(BoundaryLangevin const& dyna, BiChain const& syst, Output& output)
  {
    output.finalChainDisplay(syst.configuration(), syst.externalForce());
    if (output.doComputeCorrelations())
      output.finalDisplayAutocorrelations();
    output.displayFinalFlow(dyna.temperature(), dyna.deltaTemperature());
  }
  
  // ------------ TriChain -----------------

  
  void samplePositions(BoundaryLangevin& dyna, TriChain& syst)
  {
    cout << "Sampling the positions..." ;cout.flush();
    if (syst.isOfFixedVolum())
      for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
        syst.getParticle(iOfParticle).position(0) = 0;
    else{
      double alpha, localTemp, localBending;
      for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      {
        alpha = iOfParticle / (double) syst.nbOfParticles();
        localTemp = (1-alpha) * dyna.temperatureLeft() + alpha * dyna.temperatureRight();
        
        localBending = syst.drawPotLaw(1/localTemp);
        double position1 = (iOfParticle>0)?syst.getParticle(iOfParticle-1).position(0):0;
        double position2 = (iOfParticle>1)?syst.getParticle(iOfParticle-2).position(0):0;
        syst.getParticle(iOfParticle).position(0) = -position2 + 2 * position1 + localBending;
        syst.getParticle(iOfParticle).momentum() = syst.drawMomentum(1/localTemp, syst.getParticle(iOfParticle).mass());
      }
    }
  }  
  void writeFinalOutput(BoundaryLangevin const& dyna, TriChain const& syst, Output& output)
  {
    output.finalChainDisplay(syst.configuration(), syst.externalForce());
    if (output.doComputeCorrelations())
      output.finalDisplayAutocorrelations();
    output.displayFinalFlow(dyna.temperature(), dyna.deltaTemperature(), dyna.tauBending(), dyna.xi());
  }
  


}
