#include "Simulation.hpp"

namespace simol {
  
  //void sampleSystem(Dynamics& /*dyna*/, System& /*syst*/)
  //{throw std::invalid_argument("sampleSystem not defined");}
  
  //void writeFinalOutput(Dynamics const& dyna, System const& syst, Output& output);

  
  //###### ISOLATED #####
  
  void sampleSystem(Hamiltonian& /*dyna*/, Isolated& /*syst*/)
  {}
  
  void sampleSystem(UniformStochasticDynamics& dyna, Isolated& syst)
  {
    syst.getParticle(0).momentum() = syst.drawMomentum(dyna.beta(), syst.getParticle(0).mass());
    syst.getParticle(0).position(0) = syst.drawPotLaw(dyna.beta());
  }
  
  //template <>
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
  
  void writeFinalOutput(Hamiltonian const& dyna, Isolated const& syst, Output& output)
  {
    output.finalDisplay(syst.configuration(), dyna.externalForce());
    if (output.doComputeCorrelations())
      output.finalDisplayAutocorrelations();
    output.displayFinalVelocity(0, dyna.externalForce(0), output.velocityCV_->nbOfFourier(), output.velocityCV_->nbOfHermite());
  }
  
  //###### FLUID ######
  
  void writeFinalOutput(Dynamics const& dyna, Fluid const& syst, Output& output)
  {
    output.finalDisplay(syst.configuration(), dyna.externalForce());
    if (output.doComputeCorrelations())
      output.finalDisplayAutocorrelations();
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

  //template <>
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
  
  void writeFinalOutput(BoundaryLangevin const& dyna, TriChain const& syst, Output& output)
  {
    output.finalDisplay(syst.configuration(), dyna.externalForce());
    if (output.doComputeCorrelations())
      output.finalDisplayAutocorrelations();
    output.displayFinalFlow(dyna.temperature(), dyna.deltaTemperature(), dyna.tauBending(), dyna.xi());
  }
  
  // ---------------- Classified by Dynamics -------------------
  
  // ###### DPDE ######
  
  void sampleSystem(DPDE& dyna, Isolated& syst)
  {
    syst.getParticle(0).momentum() = syst.drawMomentum(dyna.beta(), syst.getParticle(0).mass());
    syst.getParticle(0).position(0) = 0;
    syst.getParticle(0).internalEnergy() = 1;  // TO DO : il faudra ici tirer selon la bonne loi ?
  }
  
  
  
  
  
  
  
  
}