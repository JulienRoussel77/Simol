#include "BoundaryLangevin.hpp"

namespace simol
{
  //#### BoundaryLangevin ####

  ///
  ///Constructor for Langevin dynamics on chains, where there is a thermostat at each end
  BoundaryLangevin::BoundaryLangevin(Input const& input):
    StochasticDynamics(input),
    betaLeft_(input.betaLeft()),
    betaRight_(input.betaRight()),
    temperatureLeft_(1/betaLeft_),
    temperatureRight_(1/betaRight()),
    gamma_(input.gamma()),
    tauBending_(input.tauBending())
  {
    cout << "deltaTemperature = " << deltaTemperature() << endl;
  }

  ///
  ///Read-only access for the inverse temperature at the left end
  const double& BoundaryLangevin::betaLeft() const {return betaLeft_;}
   ///
  ///Read-only access for the inverse temperature at the right end
  const double& BoundaryLangevin::betaRight() const {return betaRight_;}
  ///
  ///Returns the mean of the temperatures at the 2 ends
  double BoundaryLangevin::temperature() const
  {
    return (temperatureLeft_ + temperatureRight_) / 2;
  }
  ///
  ///Read-only access for the temperature at the left end
  const double& BoundaryLangevin::temperatureLeft() const {return temperatureLeft_;}
  ///
  ///Read-only access for the temperature at the right end
  const double& BoundaryLangevin::temperatureRight() const {return temperatureRight_;}

  ///Returns the difference between the mean temperature and the one at the left end
  ///This is equal to {eta}
  double BoundaryLangevin::deltaTemperature() const
  {
    return (temperatureLeft_ - temperatureRight_)/2;
  }
  ///
  ///Read-only accessor of the intensity of the O-U process
  const double& BoundaryLangevin::gamma() const {return gamma_;}
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
