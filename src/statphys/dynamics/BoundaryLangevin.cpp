#include "simol/statphys/dynamics/BoundaryLangevin.hpp"

namespace simol
{
  //#### BoundaryLangevin ####

  ///
  ///Constructor for Langevin dynamics on chains, where there is a thermostat at each end
  BoundaryLangevin::BoundaryLangevin(Input const& input):
    LangevinBase(input),
    deltaTemperature_(input.deltaTemperature()),
    temperatureLeft_(input.temperature() + deltaTemperature_),
    temperatureRight_(input.temperature() - deltaTemperature_),
    betaLeft_(1 / temperatureLeft_),
    betaRight_(1 / temperatureRight_),
    tauBending_(input.tauBending())
  {}

  ///
  ///Read-only access for the inverse temperature at the left end
  const double& BoundaryLangevin::betaLeft() const {return betaLeft_;}
  ///
  ///Read-only access for the inverse temperature at the right end
  const double& BoundaryLangevin::betaRight() const {return betaRight_;}
  ///
  ///Read-only access for the temperature at the left end
  const double& BoundaryLangevin::temperatureLeft() const {return temperatureLeft_;}
  ///
  ///Read-only access for the temperature at the right end
  const double& BoundaryLangevin::temperatureRight() const {return temperatureRight_;}

  ///
  ///Returns the amplitude of the brownian motion at the left end
  double BoundaryLangevin::sigmaLeft() const
  {
    return sqrt(2 * gamma_ / betaLeft_);
  }
  ///
  ///Returns the amplitude of the brownian motion at the right end
  double BoundaryLangevin::sigmaRight() const
  {
    return sqrt(2 * gamma_ / betaRight_);
  }
  ///
  ///Returns the bending constrain that is added on the right end of the chain
  const double& BoundaryLangevin::tauBending() const {return tauBending_;}
  ///
  ///Integrates the bending constraint on the (last) particle pair
  void BoundaryLangevin::bending(Particle& particle1, Particle& particle2) const
  {
    particle1.force(0) -= tauBending();
    particle2.force(0) += tauBending();
  }
  
  ///Applies the generator of this dynamics to the basis functions of the CV
  ///Evaluate at the current state of "conifguration"
  void BoundaryLangevin::computeGeneratorOnBasis(shared_ptr<CVBasis> /*cvBasis*/, System const& /*syst*/) const
  {
    /*cvBasis->generatorOnBasisValues_ = DVec::Zero(cvBasis->totalNbOfElts());
    //DVec result = DVec::Zero(nbOfFunctions());
    for (int iOfFunction = 0; iOfFunction < cvBasis->totalNbOfElts(); iOfFunction++)
    {
      //for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      cvBasis->generatorOnBasisValues_(iOfFunction) += syst.momenta().dot(cvBasis->gradientQ(syst, iOfFunction))
                              + syst.forces().dot(cvBasis->gradientP(syst, iOfFunction));
      //if(false)
      cvBasis->generatorOnBasisValues_(iOfFunction) += gamma() * (- syst(0).momentum().dot(cvBasis->basis_->gradientP(syst, iOfFunction).col(0))
                                      + cvBasis->basis_->laplacianP(syst, 0, iOfFunction) / betaLeft()
                                      - syst(syst.nbOfParticles() - 1).momentum().dot( cvBasis->basis_->gradientP(syst, iOfFunction).col(syst.nbOfParticles() - 1))
                                      + cvBasis->basis_->laplacianP(syst, syst.nbOfParticles() - 1, iOfFunction) / betaRight());
    }*/
  }
  
void BoundaryLangevin::computeProfileBiChain(Output& output, System const& syst, long int iOfStep) const
  {    
    double harOmega = syst.pairPotential().harmonicFrequency();
    double nu = syst(0).mass() * harOmega / output.constGamma_;
    
    //double alpha2 = pow(pairPotential().harmonicFrequency() / output.constGamma_, 2);
    
    static bool outbool = true;
    if (outbool)
      cout << "refFlux = " << nu * output.constDeltaTemperature_ * harOmega / (1 + pow(nu, 2)) << endl;
    outbool = false;
    
    output.obsSumFlow().currentValue() = 0;
    output.obsModiFlow().currentValue() = nu * output.constDeltaTemperature_ * harOmega / (1 + pow(nu, 2));
    //cout << nu << " " << output.constGamma_ << " " << output.constDeltaTemperature_ << " " << harOmega << endl;
    int midNb = (syst.nbOfParticles() - 1) / 2;
    output.obsMidFlow().currentValue() = syst.heatFlow(midNb);
    //cout << "---------obsModiFlow = " << output.obsModiFlow().currentValue() << endl;
    
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
    {
      double dist = 0;
      double potTempBot=1;
      double potTempTop=0;
      
      double flow = 0;
      double modiFlow = 0;
      double speedLeft, speedRight;
      
      if (iOfParticle != syst.nbOfParticles()-1)
      {
        flow = syst.heatFlow(iOfParticle);
        output.obsSumFlow().currentValue() += flow;
        
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
        
        modiFlow = -(speedRight - speedLeft) * (syst(iOfParticle).energyGrad(0) - harmonicForce) / (2 * (1+pow(nu, 2)));
        //cout << output.obsModiFlow().currentValue() << " + " << modiFlow << " = " << output.obsModiFlow().currentValue() + modiFlow << endl;
        output.obsModiFlow().currentValue() += modiFlow;
      }

      output.appendBendistProfile(dist , iOfStep, iOfParticle);
      output.appendPotTempTopProfile(potTempTop, iOfStep, iOfParticle);
      output.appendPotTempBotProfile(potTempBot, iOfStep, iOfParticle);
      output.appendKinTempProfile(2 * syst(iOfParticle).kineticEnergy(), iOfStep, iOfParticle);
      output.appendFlowProfile(flow, iOfStep, iOfParticle);
      output.appendModiFlowProfile(modiFlow, iOfStep, iOfParticle);
    }
    output.obsSumFlow().currentValue() /= (syst.nbOfParticles() - 1.);
    
  }
  
  void BoundaryLangevin::computeProfileTriChain(Output& output, System const& syst, long int iOfStep) const
  {
    output.obsSumFlow().currentValue() = 0;
    
    // modiFlow expression is easier for pair nbOfParticles
    
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
    {
      double bending = 0;
      double flow = 0;
      double potTempTop = 0;
      double potTempBot = 0;


      if (iOfParticle == 0)
      {
        bending = syst(-1).position(0) - 2 * syst(0).position(0) + syst(1).position(0);
        flow = gamma() * (parameters().temperatureLeft() - 2 * syst(0).kineticEnergy())
               - syst(-1).energyGrad(0) * syst(0).momentum(0);
        potTempTop = pow(- syst(-1).energyGrad(0) + 2 * syst(0).energyGrad(0) - syst(1).energyGrad(0), 2);
        potTempBot = syst(-1).energyLapla() + 4 * syst(0).energyLapla() + syst(1).energyLapla();
      }
      else if (iOfParticle < syst.nbOfParticles() - 1)
      {
        // bending is k_iOfParticle
        bending = syst(iOfParticle - 1).position(0) - 2 * syst(iOfParticle).position(0) + syst(iOfParticle + 1).position(0);
        // flow is j_iOfParticle
        flow = - syst(iOfParticle - 1).energyGrad(0) * syst(iOfParticle).momentum(0)
               + syst(iOfParticle).energyGrad(0) * syst(iOfParticle - 1).momentum(0);

        potTempTop = pow(- syst(iOfParticle - 1).energyGrad(0) + 2 * syst(iOfParticle).energyGrad(0) - syst(iOfParticle + 1).energyGrad(0), 2);
        potTempBot = syst(iOfParticle - 1).energyLapla() + 4 * syst(iOfParticle).energyLapla() + syst(iOfParticle + 1).energyLapla();
        output.obsSumFlow().currentValue() += flow;
      }
      else
      {
        bending = nan("");
        flow = - syst(iOfParticle - 1).energyGrad(0) * syst(iOfParticle).momentum(0);
        potTempTop = pow(- syst(iOfParticle - 1).energyGrad(0), 2);
        potTempBot = syst(iOfParticle - 1).energyLapla();
      }
      int midNb = (syst.nbOfParticles() - 1) / 2;
      if (iOfParticle == midNb)
        output.obsMidFlow().currentValue() = flow;

      output.appendBendistProfile(bending , iOfStep, iOfParticle);
      output.appendKinTempProfile(2 * syst(iOfParticle).kineticEnergy(), iOfStep, iOfParticle);
      output.appendPotTempTopProfile(potTempTop, iOfStep, iOfParticle);
      output.appendPotTempBotProfile(potTempBot, iOfStep, iOfParticle);
      output.appendFlowProfile(flow, iOfStep, iOfParticle);
    }
    output.obsSumFlow().currentValue() /= (syst.nbOfParticles() - 2.);
    //cout << output.obsSumFlow().currentValue() << endl;
  }
  
  
  
  
  //------------------------ Constrained Boundary Langevin ----------------------
  
  
  ///Constructor for the Constrained Langevin Dynamics with thermostats everywhere
  ConstrainedBoundaryLangevin::ConstrainedBoundaryLangevin(Input const& input):
    BoundaryLangevin(input),
    flux_(input.flux())
  {
    cout << "flux_ = " << flux_ << endl;
  }

}
