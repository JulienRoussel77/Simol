#include "simol/statphys/output/Output.hpp"

namespace simol
{
//   void Output::displayChainPositions(System const& syst, long int iOfStep)
//   {
//     outBeam() << iOfStep * timeStep()
//               << " " << syst(0).position() - 2 * syst(1).position() + syst(2).position()
//               //<< " " << syst(0).position() - 2*syst(1).position() + syst(2).position()
//               << " " << syst((nbOfParticles_ - 2) / 4).position() - 2 * syst((nbOfParticles_ - 2) / 4 + 1).position() + syst((nbOfParticles_ - 2) / 4 + 2).position()
//               << " " << syst((nbOfParticles_ - 2) / 2).position() - 2 * syst((nbOfParticles_ - 2) / 2 + 1).position() + syst((nbOfParticles_ - 2) / 2 + 2).position()
//               << " " << syst(3 * (nbOfParticles_ - 2) / 4).position() - 2 * syst(3 * (nbOfParticles_ - 2) / 4 + 1).position() + syst(3 * (nbOfParticles_ - 2) / 4 + 2).position()
//               << " " << syst(nbOfParticles_ - 3).position() - 2 * syst(nbOfParticles_ - 2).position() + syst(nbOfParticles_ - 1).position()
//               << endl;
//   }
// 
//   void Output::displayChainMomenta(System const& syst, long int iOfStep)
//   {
//     outChainVelocities() << iOfStep * timeStep()
//                          << " " << syst(0).momentum()
//                          << " " << syst(nbOfParticles_ / 4).momentum()
//                          << " " << syst(nbOfParticles_ / 2).momentum()
//                          << " " << syst(3 * nbOfParticles_ / 4).momentum()
//                          << " " << syst(nbOfParticles_ - 1).momentum()
//                          << endl;
//   }
  


  void Output::displayInstantProfile(System const& syst, long int iOfStep)
  {
    writeInstantProfile(syst, outInstantProfile(), iOfStep);
  }
  
  void Output::displayProfile(long int iOfStep)
  {
    writeProfile(outProfile(), iOfStep);
  }
  
  void Output::writeInstantProfile(System const& syst, ofstream & out_, long int iOfStep)
  {
    assert(out_ && out_.is_open());
    for (int i = 0; i < 6; i ++)
    {
      int iOfParticle = 1 + (i * (nbOfParticles()-2))/6;
      out_ << iOfStep * timeStep() << " "
           << iOfParticle << " "
           << bendistProfile_.lastValue(iOfParticle) << " "
           << syst(iOfParticle).momentum(0) << " "
           << fluxProfile_.lastValue(iOfParticle) << " "
           << endl;
    }
  }

  void Output::writeProfile(ofstream & out_, long int iOfStep)
  {
    assert(out_ && out_.is_open());
    for (int iOfParticle = 0; iOfParticle < nbOfParticles_; iOfParticle++)
      out_ << iOfStep * timeStep() << " "
           << iOfParticle << " "
           << bendistProfile_.mean(iOfParticle) << " "
           << bendistProfile_.asymptoticVariance(iOfParticle) / sqrt(iOfStep * timeStep()) << " "
           << fluxProfile_.mean(iOfParticle) << " "
           << fluxProfile_.asymptoticVariance(iOfParticle) / sqrt(iOfStep * timeStep()) << " "
           << kinTempProfile_.mean(iOfParticle) << " "
           << kinTempProfile_.asymptoticVariance(iOfParticle) / sqrt(iOfStep * timeStep()) << " "
           << potTempTopProfile_.mean(iOfParticle) / potTempBotProfile_.mean(iOfParticle) << " "
           << extFluxProfile_.mean(iOfParticle) << " "
           << extFluxProfile_.asymptoticVariance(iOfParticle) / sqrt(iOfStep * timeStep()) << " "
           << endl;
  }

  void Output::finalChainDisplay()
  {
    writeProfile(outFinalProfile(), nbOfSteps());
  }

  void Output::displayFinalFlux(double parameter1, double parameter2, double parameter3)
  {
    cout << "outFinalFlux_ : " <<  std::left << setw(10) << finalTime()
                     << " " << setw(5) << timeStep()
                     << " " << setw(6) << nbOfParticles()
                     << " " << setw(6) << parameters_.temperature()
                     << " " << setw(6) << parameters_.eta()
                     << " " << setw(6) << parameters_.nu()
                     << " " << setw(6) << parameter1
                     << " " << setw(6) << parameter2
                     << " " << setw(6) << parameter3
                     << " " << setw(12) << obsMidFlux_->mean()
                     << " " << setw(12) << obsMidFlux_->asymptoticVariance()
                     << " " << setw(12) << obsMidFlux_->asyvarOfAsyvar()
                     << " " << setw(12) << obsSumFlux_->mean()
                     << " " << setw(12) << obsSumFlux_->asymptoticVariance()     //#12
                     << " " << setw(12) << obsSumFlux_->asyvarOfAsyvar()
                     << " " << setw(12) << obsModiFlux_->mean()
                     << " " << setw(12) << obsModiFlux_->asymptoticVariance()
                     << " " << setw(12) << obsModiFlux_->asyvarOfAsyvar()
         << std::endl;

    //cout << "displayFinalFlux(double temperature, double delta_temperature, double tau)";
    if (doFinalFlux_)
    {      
      outFinalFlux() << std::left << setw(10) << finalTime()
                     << " " << setw(5) << timeStep()
                     << " " << setw(6) << nbOfParticles()
                     << " " << setw(4) << parameters_.temperature()
                     << " " << setw(4) << parameters_.eta()
                     << " " << setw(6) << parameters_.nu()
                     << " " << setw(6) << parameter1
                     << " " << setw(6) << parameter2
                     << " " << setw(6) << parameter3
                     << " " << setw(12) << obsMidFlux_->mean()
                     << " " << setw(12) << obsMidFlux_->asymptoticVariance()
                     << " " << setw(12) << obsMidFlux_->asyvarOfAsyvar()
                     << " " << setw(12) << obsSumFlux_->mean()
                     << " " << setw(12) << obsSumFlux_->asymptoticVariance()     //#12
                     << " " << setw(12) << obsSumFlux_->asyvarOfAsyvar()
                     << " " << setw(12) << obsModiFlux_->mean()
                     << " " << setw(12) << obsModiFlux_->asymptoticVariance()
                     << " " << setw(12) << obsModiFlux_->asyvarOfAsyvar()
                     << std::endl;
    }
  }
  
  /// For a FPU potential a,b,c : p1 = b, p2 = c and p3 = fitted stiffness
  void Output::displayFinalChainLagrangeMultiplier(double parameter1, double parameter2, double parameter3)
  {
    cout << "outFinalChainLagrangeMultiplier_ : " <<  std::left << setw(10) << finalTime()
                     << " " << setw(5) << timeStep()
                     << " " << setw(6) << nbOfParticles()
                     << " " << setw(6) << parameters_.temperature()
                     << " " << setw(6) << parameters_.flux()
                     << " " << setw(6) << parameter1
                     << " " << setw(6) << parameter2
                     << " " << setw(6) << parameter3
                     << " " << setw(12) << obsLagrangeMultiplier().mean()
                     << " " << setw(12) << obsLagrangeMultiplier().asymptoticVariance()
                     << " " << setw(12) << obsLagrangeMultiplier().asyvarOfAsyvar();

    //cout << "displayFinalFlux(double temperature, double delta_temperature, double tau)";
    if (doFinalFlux_)
    {
      outFinalLagrangeMultiplier() << std::left << setw(10) << finalTime()
                     << " " << setw(5) << timeStep()
                     << " " << setw(6) << nbOfParticles()
                     << " " << setw(4) << parameters_.temperature()
                     << " " << setw(4) << parameters_.flux()
                     << " " << setw(6) << parameter1
                     << " " << setw(6) << parameter2
                     << " " << setw(6) << parameter3
                     << " " << setw(12);
      obsLagrangeMultiplier().displayFinalValues(outFinalLagrangeMultiplier());  
      //obsMidFlux().displayFinalValues(outFinalLagrangeMultiplier());             
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

  void Output::appendFluxProfile(double value, long int iOfStep, int iOfParticle)
  {
    fluxProfile_.append(value, iOfStep, iOfParticle);
  }
  
  void Output::appendModiFluxProfile(double value, long int iOfStep, int iOfParticle)
  {
    modiFluxProfile_.append(value, iOfStep, iOfParticle);
  }
  
  void Output::appendExtFluxProfile(double value, long int iOfStep, int iOfParticle)
  {
    extFluxProfile_.append(value, iOfStep, iOfParticle);
  }

}
