#include "simol/statphys/dynamics/BoundaryLangevin.hpp"

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
  
  ///Applies the generator of this dynamics to the basis functions of the CV
  ///Evaluate at the current state of "conifguration"
  void BoundaryLangevin::computeGeneratorOnBasis(shared_ptr<CVBasis> /*cvBasis*/, System const& /*syst*/) const
  {}
  
void BoundaryLangevin::computeProfileBiChain(Output& output, System const& syst, long int iOfStep) const
  {    
    double harOmega = syst.pairPotential().harmonicFrequency();
    double nu = syst(0).mass() * harOmega / gamma();
    
    //double alpha2 = pow(pairPotential().harmonicFrequency() / output.constGamma_, 2);
    
    static bool outbool = true;
    if (outbool)
      cout << "refFlux = " << nu * deltaTemperature() * harOmega / (1 + pow(nu, 2)) << endl;
    outbool = false;
    
    output.obsSumFlux().currentValue() = 0;
    output.obsModiFlux().currentValue() = nu * deltaTemperature()* harOmega / (1 + pow(nu, 2));
    //cout << nu << " " << output.constGamma_ << " " << output.constDeltaTemperature_ << " " << harOmega << endl;
    int midNb = (syst.nbOfParticles() - 1) / 2;
    output.obsMidFlux().currentValue() = syst.heatFlux(midNb);
    //cout << "---------obsModiFlux = " << output.obsModiFlux().currentValue() << endl;
    
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
    {
      double dist = 0;
      double potTempBot=1;
      double potTempTop=0;
      
      double flux = 0;
      double modiFlux = 0;
      double speedLeft, speedRight;
      
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
        //cout << output.obsModiFlux().currentValue() << " + " << modiFlux << " = " << output.obsModiFlux().currentValue() + modiFlux << endl;
        output.obsModiFlux().currentValue() += modiFlux;
      }
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
  
  void BoundaryLangevin::computeProfileTriChain(Output& output, System const& syst, long int iOfStep) const
  {
    output.obsSumFlux().currentValue() = 0;
    
    // modiFlux expression is easier for pair nbOfParticles
    
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
    {
      double bending = 0;
      double flux = 0;
      double potTempTop = 0;
      double potTempBot = 0;


      if (iOfParticle == 0)
      {
        bending = syst(-1).position(0) - 2 * syst(0).position(0) + syst(1).position(0);
        flux = gamma() * (parameters().temperatureLeft() - 2 * syst(0).kineticEnergy())
               - syst(-1).energyGrad(0) * syst(0).momentum(0);
        potTempTop = pow(- syst(-1).energyGrad(0) + 2 * syst(0).energyGrad(0) - syst(1).energyGrad(0), 2);
        potTempBot = syst(-1).energyLapla() + 4 * syst(0).energyLapla() + syst(1).energyLapla();
      }
      else if (iOfParticle < syst.nbOfParticles() - 1)
      {
        // bending is k_iOfParticle
        bending = syst(iOfParticle - 1).position(0) - 2 * syst(iOfParticle).position(0) + syst(iOfParticle + 1).position(0);
        // flux is j_iOfParticle
        flux = - syst(iOfParticle - 1).energyGrad(0) * syst(iOfParticle).momentum(0)
               + syst(iOfParticle).energyGrad(0) * syst(iOfParticle - 1).momentum(0);

        potTempTop = pow(- syst(iOfParticle - 1).energyGrad(0) + 2 * syst(iOfParticle).energyGrad(0) - syst(iOfParticle + 1).energyGrad(0), 2);
        potTempBot = syst(iOfParticle - 1).energyLapla() + 4 * syst(iOfParticle).energyLapla() + syst(iOfParticle + 1).energyLapla();
        output.obsSumFlux().currentValue() += flux;
      }
      else
      {
        bending = nan("");
        flux = - syst(iOfParticle - 1).energyGrad(0) * syst(iOfParticle).momentum(0);
        potTempTop = pow(- syst(iOfParticle - 1).energyGrad(0), 2);
        potTempBot = syst(iOfParticle - 1).energyLapla();
      }
      int midNb = (syst.nbOfParticles() - 1) / 2;
      if (iOfParticle == midNb)
        output.obsMidFlux().currentValue() = flux;

      output.appendBendistProfile(bending , iOfStep, iOfParticle);
      output.appendKinTempProfile(2 * syst(iOfParticle).kineticEnergy(), iOfStep, iOfParticle);
      output.appendPotTempTopProfile(potTempTop, iOfStep, iOfParticle);
      output.appendPotTempBotProfile(potTempBot, iOfStep, iOfParticle);
      output.appendFluxProfile(flux, iOfStep, iOfParticle);
    }
    output.obsSumFlux().currentValue() /= (syst.nbOfParticles() - 2.);
    //cout << output.obsSumFlux().currentValue() << endl;
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
    
    //double flux = syst.leftHeatFlux();
    
    
    // Becareful here we assume that all the particles share the same mass !
    
    //syst.enforceConstraint(syst.lagrangeMultiplier(), flux(), parameters());
    
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      updateMomentum(syst(iOfParticle));
    syst.enforceConstraint(syst.lagrangeMultiplier(), flux(), parameters());
    
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      updatePosition(syst(iOfParticle));
    syst.computeAllForces();

    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      updateMomentum(syst(iOfParticle));
    syst.enforceConstraint(syst.lagrangeMultiplier(), flux(), parameters());
    
    updateOrsteinUhlenbeck(syst.getParticle(0), betaLeft(), timeStep()/2);
    updateOrsteinUhlenbeck(syst.getParticle(syst.nbOfParticles() - 1), betaRight(), timeStep()/2);
    
    // Becareful here we assume that all the particles share the same mass !
    
    //syst.enforceConstraint(syst.lagrangeMultiplier(), flux(), parameters());
    
    syst.lagrangeMultiplier() /= timeStep();
  }

}
