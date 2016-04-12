#include "BoundaryLangevin.hpp"

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

}
