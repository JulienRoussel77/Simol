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
  void BoundaryLangevin::computeGeneratorOnBasis(CVBasis& cvBasis, System const& syst) const
  {
    cvBasis.generatorOnBasisValues_ = DVec(cvBasis.totalNbOfElts(), 0);
    //Vector<double> result = Vector<double>::Zero(nbOfFunctions());
    for (int iOfFunction = 0; iOfFunction < cvBasis.totalNbOfElts(); iOfFunction++)
    {
      for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
        cvBasis.generatorOnBasisValues_(iOfFunction) += dot(syst(iOfParticle).momentum(), cvBasis.basis_->gradientQ(syst, iOfParticle, iOfFunction))
                               + dot(syst(iOfParticle).force(), cvBasis.basis_->gradientP(syst, iOfParticle, iOfFunction));
      //if(false)
      cvBasis.generatorOnBasisValues_(iOfFunction) += gamma() * (- dot(syst(0).momentum(), cvBasis.basis_->gradientP(syst, 0, iOfFunction))
                                      + cvBasis.basis_->laplacianP(syst, 0, iOfFunction) / betaLeft()
                                      - dot(syst(syst.nbOfParticles() - 1).momentum(), cvBasis.basis_->gradientP(syst, syst.nbOfParticles() - 1, iOfFunction))
                                      + cvBasis.basis_->laplacianP(syst, syst.nbOfParticles() - 1, iOfFunction) / betaRight());
    }
  }
  
  
  void BoundaryLangevin::computeProfileBiChain(Output& output, System const& syst, long int iOfStep) const
  {
    //output.obsSumFlow().currentValue() = 0;
    output.obsModiFlow().currentValue() = output.constGamma_ * output.constDeltaTemperature_ / 2;
    //cout << output.constGamma_ << " " << output.constDeltaTemperature_<< endl;
    assert(syst.nbOfParticles() % 2 == 0);
    int midNb = (syst.nbOfParticles() - 1) / 2;
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
    {
      //Particle& particle = configuration_[iOfParticle];
      double dist = 0;
      double distPrev = 0;
      double distNext = 0;
      double flow = 0;
      double potTempTop = 0;
      double potTempBot = 0;
      

      if (iOfParticle == 0)
      {
        dist = syst(0).position(0) - syst(-1).position(0);
        flow = gamma() * (temperatureLeft() - 2 * syst(0).kineticEnergy());
        //distNe = syst(1).position(0) - syst(0).position(0);
        //output.obsModiFlow().currentValue() += (syst(0).momentum(0) + distNe) / 4 * (syst(0).energyGrad(0) - dist);
        //cout << iOfParticle << " : " << syst(0).energyGrad(0) << " " << dist << endl;
      }
      else
      {
        // dist is r_iOfParticle
        dist = syst(iOfParticle).position(0) - syst(iOfParticle - 1).position(0);
        if (iOfParticle != 1)
          distPrev = syst(iOfParticle-1).position(0) - syst(iOfParticle-2).position(0);
        // flow is j_iOfParticle
        flow = - syst(iOfParticle).energyGrad(0) * syst(iOfParticle - 1).momentum(0);
        if (iOfParticle != syst.nbOfParticles() - 1)
        {
          output.obsSumFlow().currentValue() += flow;
          if (iOfParticle == midNb)
            output.obsMidFlow().currentValue() = flow;
          
          distNext = syst(iOfParticle+1).position(0) - syst(iOfParticle).position(0);
          
          if (iOfParticle == 1)
            output.obsModiFlow().currentValue() += (distNext - syst(0).momentum(0)) / 4 * (syst(iOfParticle).energyGrad(0) - dist);
          // iOfParticle != 0, 1, N-1
          else
            output.obsModiFlow().currentValue() += (distNext - distPrev) / 4 * (syst(iOfParticle).energyGrad(0) - dist);
        }
        else
          output.obsModiFlow().currentValue() -= (syst(iOfParticle).momentum(0) + distPrev) / 4 * (syst(iOfParticle).energyGrad(0) - dist);
        
        //cout << iOfParticle << " : " << syst(iOfParticle).energyGrad(0) << " " << dist << endl;
      }

      if (iOfParticle == syst.nbOfParticles() - 1)
      {
        potTempTop = pow(syst(iOfParticle).energyGrad(0), 2);
        potTempBot = syst(iOfParticle).energyLapla();
      }
      else
      {
        potTempTop = pow(syst(iOfParticle).energyGrad(0) - syst(iOfParticle + 1).energyGrad(0), 2);
        potTempBot = syst(iOfParticle).energyLapla() + syst(iOfParticle + 1).energyLapla();
      }
      
      

      output.appendBendistProfile(dist , iOfStep, iOfParticle);
      output.appendKinTempProfile(2 * syst(iOfParticle).kineticEnergy(), iOfStep, iOfParticle);
      output.appendPotTempTopProfile(potTempTop, iOfStep, iOfParticle);
      output.appendPotTempBotProfile(potTempBot, iOfStep, iOfParticle);
      output.appendFlowProfile(flow, iOfStep, iOfParticle);
    }
    output.obsSumFlow().currentValue() /= (syst.nbOfParticles() - 2.);
    //output.obsModiFlow().currentValue() /= syst.nbOfParticles();
    //cout << output.obsSumFlow().currentValue() << endl;
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
        flow = gamma() * (temperatureLeft() - 2 * syst(0).kineticEnergy())
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

}
