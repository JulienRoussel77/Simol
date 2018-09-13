#include "simol/dynamics/BoundaryLangevin.hpp"

namespace simol
{
  //#### BoundaryLangevin ####

  ///
  ///Constructor for Langevin dynamics on chains, where there is a thermostat at each end
  BoundaryLangevin::BoundaryLangevin(Input const& input):
    LangevinBase(input)
  {}

  ///
  ///Read-only access for the inverse temperature at the left end
  const double& BoundaryLangevin::betaLeft() const {return parameters_.betaLeft();}
  ///
  ///Read-only access for the inverse temperature at the right end
  const double& BoundaryLangevin::betaRight() const {return parameters_.betaRight();}
  ///
  ///Read-only access for the temperature at the left end
  const double& BoundaryLangevin::temperatureLeft() const {return parameters_.temperatureLeft();}
  ///
  ///Read-only access for the temperature at the right end
  const double& BoundaryLangevin::temperatureRight() const {return parameters_.temperatureRight();}

  ///
  ///Returns the amplitude of the brownian motion at the left end
  double BoundaryLangevin::sigmaLeft() const
  {
    return sqrt(2 * gamma() / parameters_.betaLeft());
  }
  ///
  ///Returns the amplitude of the brownian motion at the right end
  double BoundaryLangevin::sigmaRight() const
  {
    return sqrt(2 * gamma() / parameters_.betaRight());
  }
  ///
  ///Returns the bending constrain that is added on the right end of the chain
  const double& BoundaryLangevin::tauBending() const {return parameters_.tauBending();}
  ///
  ///Integrates the bending constraint on the (last) particle pair
  void BoundaryLangevin::bending(Particle& particle1, Particle& particle2) const
  {
    particle1.force(0) -= tauBending();
    particle2.force(0) += tauBending();
  }
  
  void BoundaryLangevin::thermalize(System& syst) const
  {
    for (int i = 0; i < syst.nbOfParticles(); i++)
      verletFirstPart(syst(i));

    syst.computeAllForces();

    for (int i = 0; i < syst.nbOfParticles(); i++)
      verletSecondPart(syst(i));

    for (int i = 0; i < syst.nbOfParticles(); i++)
    {
      double localTemperature = temperatureLeft() + i * deltaTemperature() / syst.nbOfParticles();
      updateOrsteinUhlenbeck(syst(0), 1 / localTemperature, timeStep());
    }
  }
  
  void BoundaryLangevin::simulate(System& syst) const
  {
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      verletFirstPart(syst(iOfParticle));

    syst.computeAllForces();

    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      verletSecondPart(syst(iOfParticle));

    updateOrsteinUhlenbeck(syst.getParticle(0), betaLeft(), timeStep());
    updateOrsteinUhlenbeck(syst.getParticle(syst.nbOfParticles() - 1), betaRight(), timeStep());

    if (doMomentaExchange())
      for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles() - 1; iOfParticle++)
        updateMomentaExchange(syst(iOfParticle), syst(iOfParticle + 1));
  }
  
void BoundaryLangevin::computeProfileChain(Output& output, System const& syst, long int iOfStep) const
  {    
    double harOmega = syst.pairPotential().harmonicFrequency();
    double nu = syst(0).mass() * harOmega / gamma();
    
    output.obsSumFlux().currentValue() = 0;
    output.obsModiFlux().currentValue() = nu * deltaTemperature()* harOmega / (1 + pow(nu, 2));
    
    int midNb = (syst.nbOfParticles() - 1) / 2;
    output.obsMidFlux().currentValue() = syst.heatFlux(midNb);
        
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
    {
      double dist = 0;
      double potTempBot=1;
      double potTempTop=0;
      
      double flux = 0;
      double modiFlux = 0;
      double speedLeft, speedRight;
      
      // computes the different fluxes, plus the configurational temperature
      if (iOfParticle != syst.nbOfParticles()-1)
      {
        if (iOfParticle != 0)
          flux = syst.heatFluxOnAtom(iOfParticle);
        output.obsSumFlux().currentValue() += flux;
        
        dist = syst(iOfParticle+1).position(0) - syst(iOfParticle).position(0);        // dist is r_iOfParticle
        
        potTempTop = pow(syst(iOfParticle).energyGrad(0), 2);
        potTempBot = syst(iOfParticle).energyLapla();
        double harmonicForce = syst.pairPotential().harmonicForce(dist);
        
        if (iOfParticle == 0)
          speedLeft = -syst(0).momentum(0) / syst(0).mass();
        else
          speedLeft = -nu * harOmega * (syst(iOfParticle).position(0) - syst(iOfParticle-1).position(0) - syst.pairPotential().harmonicEquilibrium());
          
        if (iOfParticle == syst.nbOfParticles()-2)
          speedRight = syst(syst.nbOfParticles()-1).momentum(0) / syst(0).mass();
        else
          speedRight = -nu * harOmega * (syst(iOfParticle+2).position(0) - syst(iOfParticle+1).position(0) - syst.pairPotential().harmonicEquilibrium());
        
        modiFlux = -(speedRight - speedLeft) * (syst(iOfParticle).energyGrad(0) - harmonicForce) / (2 * (1+pow(nu, 2)));
        
        output.obsModiFlux().currentValue() += modiFlux;
      }
      
      // updates the statstics on the the profile of some quantities along the chain
      output.appendBendistProfile(dist , iOfStep, iOfParticle);
      output.appendPotTempTopProfile(potTempTop, iOfStep, iOfParticle);
      output.appendPotTempBotProfile(potTempBot, iOfStep, iOfParticle);
      output.appendFluxProfile(flux, iOfStep, iOfParticle);
      output.appendModiFluxProfile(modiFlux, iOfStep, iOfParticle);
      output.appendKinTempProfile(2 * syst(iOfParticle).kineticEnergy(), iOfStep, iOfParticle); 
      output.appendExtFluxProfile(flux * syst.lagrangeMultiplier(), iOfStep, iOfParticle);
    }
    output.obsSumFlux().currentValue() /= (syst.nbOfParticles() - 2.);
    
  }
  
  
  
  //------------------------ Constrained Boundary Langevin ----------------------
  
  
  ///Constructor for the Constrained Langevin Dynamics with thermostats everywhere
  ConstrainedBoundaryLangevin::ConstrainedBoundaryLangevin(Input const& input):
    BoundaryLangevin(input)
  {}
  
  void ConstrainedBoundaryLangevin::simulate(System& syst) const
  {
    syst.lagrangeMultiplier() = 0;
    updateOrsteinUhlenbeck(syst.getParticle(0), betaLeft(), timeStep()/2);
    updateOrsteinUhlenbeck(syst.getParticle(syst.nbOfParticles() - 1), betaRight(), timeStep()/2);
        
    // Becareful here we assume that all the particles share the same mass !
    
    // The flux does not depend on p_1 and p_N so we don't reproject
    //syst.enforceConstraint(syst.lagrangeMultiplier(), flux(), parameters());
    
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      updateMomentum(syst(iOfParticle));
    syst.enforceConstraint(flux(), parameters(), true);
    
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      updatePosition(syst(iOfParticle));
    syst.computeAllForces();

    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      updateMomentum(syst(iOfParticle));
    syst.enforceConstraint(flux(), parameters(), true);
    
    updateOrsteinUhlenbeck(syst.getParticle(0), betaLeft(), timeStep()/2);
    updateOrsteinUhlenbeck(syst.getParticle(syst.nbOfParticles() - 1), betaRight(), timeStep()/2);
    
    // Becareful here we assume that all the particles share the same mass !
    
    //syst.enforceConstraint(syst.lagrangeMultiplier(), flux(), parameters());
    
    syst.lagrangeMultiplier() /= timeStep();
  }

}
