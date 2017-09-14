#include "simol/statphys/simulation/Simulation.hpp"

namespace simol
{




  void thermalize(BoundaryLangevin& dyna, System& syst)
  {
    //for (auto&& particle : configuration_)
    //dyna.updateBefore(particle);
    for (int i = 0; i < syst.nbOfParticles(); i++)
      dyna.verletFirstPart(syst(i));

    syst.computeAllForces();

    for (int i = 0; i < syst.nbOfParticles(); i++)
      dyna.verletSecondPart(syst(i));

    for (int i = 0; i < syst.nbOfParticles(); i++)
    {
      double localTemperature = dyna.temperatureLeft() + i * dyna.deltaTemperature() / syst.nbOfParticles();
      dyna.updateOrsteinUhlenbeck(syst(0), 1 / localTemperature, dyna.timeStep());
    }
  }
  
    //------------ BoundaryLangevin ---------------
  
  void simulate(BoundaryLangevin& dyna, System& syst)
  {
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      dyna.verletFirstPart(syst(iOfParticle));

    syst.computeAllForces();

    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      dyna.verletSecondPart(syst(iOfParticle));

    dyna.updateOrsteinUhlenbeck(syst.getParticle(0), dyna.betaLeft(), dyna.timeStep());
    dyna.updateOrsteinUhlenbeck(syst.getParticle(syst.nbOfParticles() - 1), dyna.betaRight(), dyna.timeStep());

    if (dyna.doMomentaExchange())
      for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles() - 1; iOfParticle++)
        dyna.updateMomentaExchange(syst(iOfParticle), syst(iOfParticle + 1));
  }
  
  
  void simulate(ConstrainedBoundaryLangevin& dyna, System& syst)
  {
    dyna.lagrangeMultiplier() = 0;
    //double firstMomentum = syst.getParticle(0).momentum(0);
    dyna.updateOrsteinUhlenbeck(syst.getParticle(0), dyna.betaLeft(), dyna.timeStep()/2);
    //double lastMomentum = syst.getParticle(syst.nbOfParticles()-1).momentum(0);
    dyna.updateOrsteinUhlenbeck(syst.getParticle(syst.nbOfParticles() - 1), dyna.betaRight(), dyna.timeStep()/2);
    
    //double flux = syst.leftHeatFlow();
    
    
    // Becareful here we assume that all the particles share the same mass !
    // We analiticaly determine the non-martingale part of the Lagrange multiplier
    //double alpha = exp(- dyna.gamma() / syst(0).mass() * dyna.timeStep()/2);
    //dyna.lagrangeMultiplier() += (1-alpha) * dyna.drift();
    //double trash=0;
    //syst.enforceConstraint(trash, dyna.drift());
    syst.enforceConstraint(dyna.lagrangeMultiplier(), dyna.flux(), dyna.parameters());
    
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      dyna.updateMomentum(syst(iOfParticle));
    syst.enforceConstraint(dyna.lagrangeMultiplier(), dyna.flux(), dyna.parameters());
    
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      dyna.updatePosition(syst(iOfParticle));
    syst.computeAllForces();

    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      dyna.verletSecondPart(syst(iOfParticle));
    syst.enforceConstraint(dyna.lagrangeMultiplier(), dyna.flux(), dyna.parameters());
    
    dyna.updateOrsteinUhlenbeck(syst.getParticle(0), dyna.betaLeft(), dyna.timeStep()/2);
    dyna.updateOrsteinUhlenbeck(syst.getParticle(syst.nbOfParticles() - 1), dyna.betaRight(), dyna.timeStep()/2);
    
    // Becareful here we assume that all the particles share the same mass !
    // We analiticaly retermine the non-martingale part of the Lagrange multiplier
    //dyna.lagrangeMultiplier() += (1-alpha) * dyna.drift();
    //syst.enforceConstraint(trash, dyna.drift());
    syst.enforceConstraint(dyna.lagrangeMultiplier(), dyna.flux(), dyna.parameters());
    
    dyna.lagrangeMultiplier() /= dyna.timeStep();
  }
  
  
  
  




}
