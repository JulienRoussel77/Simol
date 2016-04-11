#include "Output.hpp"

namespace simol{
  
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
  {return *outSumFlowCV_;}
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
  
  ControlVariate& Output::midFlowCV()
  {return *midFlowCV_;}
  
  ControlVariate& Output::sumFlowCV()
  {return *sumFlowCV_;}
  
  void Output::displayChainPositions(vector<Particle> const& configuration, int iOfIteration)
  {
    outBeam() << iOfIteration * timeStep() 
        << " " << configuration[0].position() - 2*configuration[1].position() + configuration[2].position()
      //<< " " << configuration[0].position() - 2*configuration[1].position() + configuration[2].position()
        << " " << configuration[(nbOfParticles_-2)/4].position() - 2*configuration[(nbOfParticles_-2)/4+1].position() + configuration[(nbOfParticles_-2)/4+2].position()
        << " " << configuration[(nbOfParticles_-2)/2].position() - 2*configuration[(nbOfParticles_-2)/2+1].position() + configuration[(nbOfParticles_-2)/2+2].position()
        << " " << configuration[3*(nbOfParticles_-2)/4].position() - 2*configuration[3*(nbOfParticles_-2)/4+1].position() + configuration[3*(nbOfParticles_-2)/4+2].position()
        << " " << configuration[nbOfParticles_-3].position() - 2*configuration[nbOfParticles_-2].position() + configuration[nbOfParticles_-1].position()
        << endl;   
  }
  
  void Output::displayChainMomenta(vector<Particle> const& configuration, int iOfIteration)
  {
    outChainVelocities() << iOfIteration * timeStep()
      << " " << configuration[0].momentum()
      << " " << configuration[nbOfParticles_/4].momentum()
      << " " << configuration[nbOfParticles_/2].momentum()
      << " " << configuration[3*nbOfParticles_/4].momentum()
      << " " << configuration[nbOfParticles_-1].momentum()
      << endl;
  }
  
   void Output::displayProfile(int iOfIteration)
  {
    writeProfile(outProfile(), iOfIteration);
  }

    void Output::writeProfile(ofstream & out_, int iOfIteration)
  {
    assert(out_.is_open());
    for (int iOfParticle = 0; iOfParticle < nbOfParticles_; iOfParticle++)
      out_ << iOfIteration * timeStep() << " " 
           << iOfParticle << " " 
          << bendistProfile_.mean(iOfParticle) << " "
          << bendistProfile_.standardDeviation(iOfParticle) / sqrt(iOfIteration * timeStep()) << " "
          << flowProfile_.mean(iOfParticle) << " "
          << flowProfile_.standardDeviation(iOfParticle) / sqrt(iOfIteration * timeStep()) << " "
          << kinTempProfile_.mean(iOfParticle) << " "
          << kinTempProfile_.standardDeviation(iOfParticle) / sqrt(iOfIteration * timeStep()) << " "
          << potTempTopProfile_.mean(iOfParticle) / potTempBotProfile_.mean(iOfParticle) << " "
          << endl;
  }
  
  void Output::finalChainDisplay(vector<Particle> const& /*configuration*/, Vector<double> const& /*externalForce*/)
  {    
    writeProfile(outFinalProfile(), nbOfIterations());
    
    midFlowCV_->postTreat(outMidFlowPT(), timeStep());
    sumFlowCV_->postTreat(outSumFlowPT(), timeStep());
  }
  
    void Output::displayFinalFlow(double temperature, double delta_temperature, double tau, double xi)
  {
    cout << "outFinalFlow_ : " <<  std::left << setw(10) << finalTime()
      << " " << setw(5) << timeStep()
      << " " << setw(6) << nbOfParticles()
      << " " << setw(4) << temperature
      << " " << setw(4) << delta_temperature
      << " " << setw(4) << tau
      << " " << setw(6) << xi
      << " " << setw(12) << midFlowCV_->meanObservable()
      << " " << setw(12) << midFlowCV_->stdDeviationObservable()
      << " " << setw(12) << sumFlowCV_->meanObservable()
      << " " << setw(12) << sumFlowCV_->stdDeviationObservable()
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
      << " " << setw(4) << tau
      << " " << setw(6) << xi
      << " " << setw(12) << midFlowCV_->meanObservable()
      << " " << setw(12) << midFlowCV_->stdDeviationObservable()
      << " " << setw(12) << sumFlowCV_->meanObservable()
      << " " << setw(12) << sumFlowCV_->stdDeviationObservable()
      << std::endl;
    }
  }
  
  void Output::appendKinTempProfile(double value, int iOfIteration, int iOfParticle)
  {
    kinTempProfile_.append(value, iOfIteration, iOfParticle);
  }
  
  void Output::appendPotTempTopProfile(double value, int iOfIteration, int iOfParticle)
  {
    potTempTopProfile_.append(value, iOfIteration, iOfParticle);
  }
  
      void Output::appendPotTempBotProfile(double value, int iOfIteration, int iOfParticle)
  {
    potTempBotProfile_.append(value, iOfIteration, iOfParticle);
  }
  
  void Output::appendBendistProfile(double value, int iOfIteration, int iOfParticle)
  {
    bendistProfile_.append(value, iOfIteration, iOfParticle);
  }
  
  void Output::appendFlowProfile(double value, int iOfIteration, int iOfParticle)
  {
    flowProfile_.append(value, iOfIteration, iOfParticle);
  }
  
}
