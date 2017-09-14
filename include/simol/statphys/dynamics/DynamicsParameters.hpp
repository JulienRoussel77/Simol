#ifndef SIMOL_DYNAMICSPARAMETERS_HPP
#define SIMOL_DYNAMICSPARAMETERS_HPP

#include "simol/statphys/input/Input.hpp"

namespace simol
{
  class DynamicsParameters
  {
  public:
    DynamicsParameters(Input const& input):
      beta_(input.beta()),
      temperature_(1/beta_),
      deltaTemperature_(input.deltaTemperature()),
      temperatureLeft_(temperature_ + deltaTemperature_),
      temperatureRight_(temperature_ - deltaTemperature_),
      betaLeft_(1/temperatureLeft_),
      betaRight_(1/temperatureRight_),
      gamma_(input.gamma()),
      tauBending_(input.tauBending()),
      interactionRatio_(input.interactionRatio()),
      bulkDriving_(input.bulkDriving())
    {}
    const double& temperature() const {return temperature_;}
    const double& temperatureLeft() const {return temperatureLeft_;}
    const double& temperatureRight() const {return temperatureRight_;}
    const double& deltaTemperature() const {return deltaTemperature_;}
    const double& beta() const {return beta_;}
    const double& betaLeft() const {return betaLeft_;}
    const double& betaRight() const {return betaRight_;}
    const double& gamma() const {return gamma_;}
    const double& tauBending() const {return tauBending_;}
    const double& interactionRatio() const {return interactionRatio_;}
    const double& bulkDriving() const {return bulkDriving_;}

  private:
    double beta_;
    double temperature_;
    double deltaTemperature_;
    double temperatureLeft_;
    double temperatureRight_;
    double betaLeft_;
    double betaRight_;
    double gamma_;
    double tauBending_;
    double interactionRatio_;
    double bulkDriving_;

  };

}

#endif