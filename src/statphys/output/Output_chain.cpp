#include "simol/statphys/output/Output.hpp"

namespace simol
{
  void Output::displayChainPositions(System const& syst, long int iOfStep)
  {
    outBeam() << iOfStep * timeStep()
              << " " << syst(0).position() - 2 * syst(1).position() + syst(2).position()
              //<< " " << syst(0).position() - 2*syst(1).position() + syst(2).position()
              << " " << syst((nbOfParticles_ - 2) / 4).position() - 2 * syst((nbOfParticles_ - 2) / 4 + 1).position() + syst((nbOfParticles_ - 2) / 4 + 2).position()
              << " " << syst((nbOfParticles_ - 2) / 2).position() - 2 * syst((nbOfParticles_ - 2) / 2 + 1).position() + syst((nbOfParticles_ - 2) / 2 + 2).position()
              << " " << syst(3 * (nbOfParticles_ - 2) / 4).position() - 2 * syst(3 * (nbOfParticles_ - 2) / 4 + 1).position() + syst(3 * (nbOfParticles_ - 2) / 4 + 2).position()
              << " " << syst(nbOfParticles_ - 3).position() - 2 * syst(nbOfParticles_ - 2).position() + syst(nbOfParticles_ - 1).position()
              << endl;
  }

  void Output::displayChainMomenta(System const& syst, long int iOfStep)
  {
    outChainVelocities() << iOfStep * timeStep()
                         << " " << syst(0).momentum()
                         << " " << syst(nbOfParticles_ / 4).momentum()
                         << " " << syst(nbOfParticles_ / 2).momentum()
                         << " " << syst(3 * nbOfParticles_ / 4).momentum()
                         << " " << syst(nbOfParticles_ - 1).momentum()
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

  void Output::finalChainDisplay()
  {
    writeProfile(outFinalProfile(), nbOfSteps());
  }

  void Output::displayFinalFlow(double parameter1, double parameter2, double parameter3)
  {
    cout << "outFinalFlow_ : " <<  std::left << setw(10) << finalTime()
                     << " " << setw(5) << timeStep()
                     << " " << setw(6) << nbOfParticles()
                     << " " << setw(6) << constTemperature_
                     << " " << setw(6) << constDeltaTemperature_
                     << " " << setw(6) << constBulkDriving_
                     << " " << setw(6) << parameter1
                     << " " << setw(6) << parameter2
                     << " " << setw(6) << parameter3
                     << " " << setw(12) << obsMidFlow_->mean()
                     << " " << setw(12) << obsMidFlow_->variance()
                     << " " << setw(12) << obsMidFlow_->varOfVar()
                     << " " << setw(12) << obsSumFlow_->mean()
                     << " " << setw(12) << obsSumFlow_->variance()     //#12
                     << " " << setw(12) << obsSumFlow_->varOfVar()
                     << " " << setw(12) << obsModiFlow_->mean()
                     << " " << setw(12) << obsModiFlow_->variance()
                     << " " << setw(12) << obsModiFlow_->varOfVar()
         << std::endl;

    //cout << "displayFinalFlow(double temperature, double delta_temperature, double tau)";
    if (doFinalFlow_)
    {      
      outFinalFlow() << std::left << setw(10) << finalTime()
                     << " " << setw(5) << timeStep()
                     << " " << setw(6) << nbOfParticles()
                     << " " << setw(4) << constTemperature_
                     << " " << setw(4) << constDeltaTemperature_
                     << " " << setw(6) << constBulkDriving_
                     << " " << setw(6) << parameter1
                     << " " << setw(6) << parameter2
                     << " " << setw(6) << parameter3
                     << " " << setw(12) << obsMidFlow_->mean()
                     << " " << setw(12) << obsMidFlow_->variance()
                     << " " << setw(12) << obsMidFlow_->varOfVar()
                     << " " << setw(12) << obsSumFlow_->mean()
                     << " " << setw(12) << obsSumFlow_->variance()     //#12
                     << " " << setw(12) << obsSumFlow_->varOfVar()
                     << " " << setw(12) << obsModiFlow_->mean()
                     << " " << setw(12) << obsModiFlow_->variance()
                     << " " << setw(12) << obsModiFlow_->varOfVar()
                     << std::endl;
    }
  }
  
  void Output::displayFinalChainLagrangeMultiplier(double parameter1, double parameter2, double parameter3)
  {
    cout << "outFinalFlow_ : " <<  std::left << setw(10) << finalTime()
                     << " " << setw(5) << timeStep()
                     << " " << setw(6) << nbOfParticles()
                     << " " << setw(6) << constTemperature_
                     << " " << setw(6) << constDeltaTemperature_
                     << " " << setw(6) << constBulkDriving_
                     << " " << setw(6) << parameter1
                     << " " << setw(6) << parameter2
                     << " " << setw(6) << parameter3
                     << " " << setw(12);
    obsLagrangeMultiplier().displayFinalValues(outFinalLagrangeMultiplier());  
    obsMidFlow().displayFinalValues(outFinalLagrangeMultiplier());

    //cout << "displayFinalFlow(double temperature, double delta_temperature, double tau)";
    if (doFinalFlow_)
    {      
      outFinalFlow() << std::left << setw(10) << finalTime()
                     << " " << setw(5) << timeStep()
                     << " " << setw(6) << nbOfParticles()
                     << " " << setw(4) << constTemperature_
                     << " " << setw(4) << constDeltaTemperature_
                     << " " << setw(6) << constBulkDriving_
                     << " " << setw(6) << parameter1
                     << " " << setw(6) << parameter2
                     << " " << setw(6) << parameter3
                     << " " << setw(12);
      obsLagrangeMultiplier().displayFinalValues(outFinalLagrangeMultiplier());  
      obsMidFlow().displayFinalValues(outFinalLagrangeMultiplier());             
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
  
  void Output::appendModiFlowProfile(double value, long int iOfStep, int iOfParticle)
  {
    modiFlowProfile_.append(value, iOfStep, iOfParticle);
  }

}
