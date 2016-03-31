#include "UniformStochasticDynamics.hpp"

namespace simol
{

  ///
  ///Constructs a purely virtual Dynamics for Dynamics with a global temperature
  UniformStochasticDynamics::UniformStochasticDynamics(Input const& input):
    StochasticDynamics(input),
    beta_(input.beta()),
    temperature_(1/beta_)
  {}


  ///
  ///Returns the temperature
  double UniformStochasticDynamics::temperature() const {return temperature_;}
  ///
  ///Read-only access for the temperature
  const double& UniformStochasticDynamics::temperatureLeft() const {return temperature_;}
  ///
  ///Read-only access for the temperature
  const double& UniformStochasticDynamics::temperatureRight() const {return temperature_;}
  ///
  ///Read-only access for the inverse temperature
  const double& UniformStochasticDynamics::beta() const {return beta_;}
  ///
  ///Read-only access for the inverse temperature
  const double& UniformStochasticDynamics::betaLeft() const {return beta_;}
  ///
  ///Read-only access for the inverse temperature
  const double& UniformStochasticDynamics::betaRight() const {return beta_;}

}
