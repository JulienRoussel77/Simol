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
      bulkDriving_(input.bulkDriving()),
      mass_(input.mass()),
      nonEqAmplitude_(input.nonEqAmplitude()),
      nbOfQModes_(input.nbOfQModes()),
      nbOfPModes_(input.nbOfPModes()),
      drift_(input.drift()),
      flux_(input.flux()),
      xi_(input.xi()),
      coupling_(input.coupling()),
      eta_(input.eta()),
      nu_(input.nu())
      //isMollified_(input.isMollified())
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
    const double& mass() const {return mass_;}
    const double& nonEqAmplitude() const {return nonEqAmplitude_;}
    const int& nbOfQModes() const {return nbOfQModes_;}
    const int& nbOfPModes() const {return nbOfPModes_;}
    const double& drift() const {return drift_;}
    const double& flux() const {return flux_;}
    const double& xi() const {return xi_;}
    const double& coupling() const {return coupling_;}
    const double& eta() const {return eta_;}
    const double& nu() const {return nu_;}
    const bool& isMollified() const {return isMollified_;}
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
    double mass_;
    double nonEqAmplitude_;
    int nbOfQModes_, nbOfPModes_;
    double drift_;
    double flux_;
    double xi_;
    double coupling_;
    double eta_;
    double nu_;
    bool isMollified_;
  };

}

#endif