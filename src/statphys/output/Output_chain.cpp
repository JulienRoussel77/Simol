#include "simol/statphys/output/Output.hpp"

namespace simol
{

  ofstream & Output::outChainVelocities()
  {return *outChainVelocities_;}
  ofstream & Output::outBeam()
  {return *outBeam_;}
  ofstream & Output::outFinalFlow()
  {return *outFinalFlow_;}
  ofstream & Output::outMidFlowCV()
  {return *outMidFlowCV_;}
  ofstream & Output::outMidFlowPT()
  {return *outMidFlowPT_;}
  ofstream & Output::outSumFlowCV()
  {return *outSumFlowCV_;}
  ofstream & Output::outSumFlowPT()
  {return *outSumFlowPT_;}
  ofstream & Output::outModiFlowCV()
  {return *outModiFlowCV_;}
  ofstream & Output::outProfile()
  {return *outProfile_;}
  ofstream & Output::outFinalProfile()
  {return *outFinalProfile_;}

  const double& Output::energyMidFlow() const
  {return energyMidFlow_;}

  double& Output::energyMidFlow()
  {return energyMidFlow_;}

  const double& Output::energySumFlow() const
  {return energySumFlow_;}

  double& Output::energySumFlow()
  {return energySumFlow_;}
  
  const double& Output::energyModiFlow() const
  {return energyModiFlow_;}

  double& Output::energyModiFlow()
  {return energyModiFlow_;}

  ControlVariate& Output::midFlowCV()
  {return *midFlowCV_;}

  ControlVariate& Output::sumFlowCV()
  {return *sumFlowCV_;}
  
  ControlVariate& Output::modiFlowCV()
  {return *modiFlowCV_;}

  void Output::displayChainPositions(vector<Particle> const& configuration, long int iOfStep)
  {
    outBeam() << iOfStep * timeStep()
              << " " << configuration[0].position() - 2 * configuration[1].position() + configuration[2].position()
              //<< " " << configuration[0].position() - 2*configuration[1].position() + configuration[2].position()
              << " " << configuration[(nbOfParticles_ - 2) / 4].position() - 2 * configuration[(nbOfParticles_ - 2) / 4 + 1].position() + configuration[(nbOfParticles_ - 2) / 4 + 2].position()
              << " " << configuration[(nbOfParticles_ - 2) / 2].position() - 2 * configuration[(nbOfParticles_ - 2) / 2 + 1].position() + configuration[(nbOfParticles_ - 2) / 2 + 2].position()
              << " " << configuration[3 * (nbOfParticles_ - 2) / 4].position() - 2 * configuration[3 * (nbOfParticles_ - 2) / 4 + 1].position() + configuration[3 * (nbOfParticles_ - 2) / 4 + 2].position()
              << " " << configuration[nbOfParticles_ - 3].position() - 2 * configuration[nbOfParticles_ - 2].position() + configuration[nbOfParticles_ - 1].position()
              << endl;
  }

  void Output::displayChainMomenta(vector<Particle> const& configuration, long int iOfStep)
  {
    outChainVelocities() << iOfStep * timeStep()
                         << " " << configuration[0].momentum()
                         << " " << configuration[nbOfParticles_ / 4].momentum()
                         << " " << configuration[nbOfParticles_ / 2].momentum()
                         << " " << configuration[3 * nbOfParticles_ / 4].momentum()
                         << " " << configuration[nbOfParticles_ - 1].momentum()
                         << endl;
  }

  void Output::displayProfile(long int iOfStep)
  {
    writeProfile(outProfile(), iOfStep);
  }

  void Output::writeProfile(ofstream & out_, long int iOfStep)
  {
    assert(out_ && out_.is_open());
    for (int iOfParticle = 0; iOfParticle < nbOfParticles_; iOfParticle++)
      out_ << iOfStep * timeStep() << " "
           << iOfParticle << " "
           << bendistProfile_.mean(iOfParticle) << " "
           << bendistProfile_.variance(iOfParticle) / sqrt(iOfStep * timeStep()) << " "
           << flowProfile_.mean(iOfParticle) << " "
           << flowProfile_.variance(iOfParticle) / sqrt(iOfStep * timeStep()) << " "
           << kinTempProfile_.mean(iOfParticle) << " "
           << kinTempProfile_.variance(iOfParticle) / sqrt(iOfStep * timeStep()) << " "
           << potTempTopProfile_.mean(iOfParticle) / potTempBotProfile_.mean(iOfParticle) << " "
           << endl;
  }

  void Output::finalChainDisplay(vector<Particle> const& /*configuration*/, Vector<double> const& /*externalForce*/)
  {
    writeProfile(outFinalProfile(), nbOfSteps());
    midFlowCV_->postTreat(outMidFlowPT(), timeStep());
    sumFlowCV_->postTreat(outSumFlowPT(), timeStep());
  }

  void Output::displayFinalFlow(double temperature, double delta_temperature, double parameter1, double parameter2)
  {
    cout << "outFinalFlow_ : " <<  std::left << setw(10) << finalTime()
                     << " " << setw(5) << timeStep()
                     << " " << setw(6) << nbOfParticles()
                     << " " << setw(4) << temperature
                     << " " << setw(4) << delta_temperature
                     << " " << setw(4) << parameter1
                     << " " << setw(6) << parameter2
                     << " " << setw(12) << midFlowCV_->meanObservable()
                     << " " << setw(12) << midFlowCV_->varObservable()
                     << " " << setw(12) << midFlowCV_->varOfVarObservable()
                     << " " << setw(12) << sumFlowCV_->meanObservable()
                     << " " << setw(12) << sumFlowCV_->varObservable()     //#12
                     << " " << setw(12) << sumFlowCV_->varOfVarObservable()
                     << " " << setw(12) << modiFlowCV_->meanObservable()
                     << " " << setw(12) << modiFlowCV_->varObservable()
                     << " " << setw(12) << modiFlowCV_->varOfVarObservable()
         << std::endl;

    //cout << "displayFinalFlow(double temperature, double delta_temperature, double tau)";
    if (doFinalFlow_)
    {
      assert(outFinalFlow().is_open());
      outFinalFlow() << std::left << setw(10) << finalTime()
                     << " " << setw(5) << timeStep()
                     << " " << setw(6) << nbOfParticles()
                     << " " << setw(4) << temperature
                     << " " << setw(4) << delta_temperature
                     << " " << setw(4) << parameter1
                     << " " << setw(6) << parameter2
                     << " " << setw(12) << midFlowCV_->meanObservable()
                     << " " << setw(12) << midFlowCV_->varObservable()
                     << " " << setw(12) << midFlowCV_->varOfVarObservable()
                     << " " << setw(12) << sumFlowCV_->meanObservable()
                     << " " << setw(12) << sumFlowCV_->varObservable()     //#12
                     << " " << setw(12) << sumFlowCV_->varOfVarObservable()
                     << " " << setw(12) << modiFlowCV_->meanObservable()
                     << " " << setw(12) << modiFlowCV_->varObservable()
                     << " " << setw(12) << modiFlowCV_->varOfVarObservable()
                     << std::endl;
    }
  }

  void Output::appendKinTempProfile(double value, long int iOfStep, int iOfParticle)
  {
    kinTempProfile_.append(value, iOfStep, iOfParticle);
  }

  void Output::appendPotTempTopProfile(double value, long int iOfStep, int iOfParticle)
  {
    potTempTopProfile_.append(value, iOfStep, iOfParticle);
  }

  void Output::appendPotTempBotProfile(double value, long int iOfStep, int iOfParticle)
  {
    potTempBotProfile_.append(value, iOfStep, iOfParticle);
  }

  void Output::appendBendistProfile(double value, long int iOfStep, int iOfParticle)
  {
    bendistProfile_.append(value, iOfStep, iOfParticle);
  }

  void Output::appendFlowProfile(double value, long int iOfStep, int iOfParticle)
  {
    flowProfile_.append(value, iOfStep, iOfParticle);
  }

}
