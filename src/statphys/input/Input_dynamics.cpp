#include "simol/statphys/input/Input.hpp"

using std::cout;
using std::endl;
using std::string;

const double defaultSeed = 0;
const double defaultGamma = 1;
const double defaultXi = 0;
const double defaultDeltaTemperature = 0;
const double defaultTauBending = 0;
const double defaultKappa = 0; 
const int    defaultMTSfrequency = 1;
const double defaultEta = 0;
const double defaultDrift= 0;
const double defaultBulkDriving = 0;
const double defaultNu = 0;
const double defaultFlux = 0;

namespace simol
{

  string Input::dynamicsName() const {return data["Dynamics"]["Name"].as<string>();}
  
  bool Input::isConstrained() const
  {
    if (data["Dynamics"]["IsConstrained"])
      return isYes(data["Dynamics"]["IsConstrained"].as<string>());
    return false;
  }
  
//   bool Input::isMollified() const
//   {
//     if (data["Dynamics"]["IsMollified"])
//       return isYes(data["Dynamics"]["IsMollified"].as<string>());
//     return false;
//   }

  double Input::gamma() const
  {
    if (data["Dynamics"]["Gamma"])
      return data["Dynamics"]["Gamma"].as<double>();
    //else throw std::runtime_error("No Gamma in the input file !");
    else return defaultGamma;
  }

  double Input::temperature() const
  {
    if (data["Dynamics"]["Temperature"])
      return data["Dynamics"]["Temperature"].as<double>();
    else if (data["Dynamics"]["Beta"])
      return 1 / data["Dynamics"]["Beta"].as<double>();
    else if (data["Dynamics"]["TemperatureLeft"] && data["Dynamics"]["TemperatureRight"])
      return (data["Dynamics"]["TemperatureLeft"].as<double>() + data["Dynamics"]["TemperatureRight"].as<double>()) / 2;
    else if (data["Dynamics"]["BetaLeft"] && data["Dynamics"]["BetaRight"])
      return .5 / data["Dynamics"]["BetaLeft"].as<double>() + .5 / data["Dynamics"]["BetaRight"].as<double>();
    else throw std::runtime_error("No given temperature!");
  }

  double Input::beta() const
  {
    return 1 / temperature();
  }

  bool Input::doMetropolis() const
  {
    if (data["Dynamics"]["Metropolis"] && data["Dynamics"]["Metropolis"].as<string>() == "yes")
      return true;
    else
      return false;
  }

  bool Input::doProjectionDPDE() const
  {
    if (data["Dynamics"]["Projection DPDE"] && data["Dynamics"]["Projection DPDE"].as<string>() == "yes")
      return true;
    else
      return false;
  }

  double Input::deltaTemperature() const
  {
    return eta() / (2*beta());
    /*if (data["Dynamics"]["DeltaTemperature"])
      return data["Dynamics"]["DeltaTemperature"].as<double>();
    else return defaultDeltaTemperature;*/
  }

  double Input::tauBending() const
  {
    if (data["Dynamics"]["Tau"])
      return data["Dynamics"]["Tau"].as<double>();
    else
      return defaultTauBending;
  }

  double Input::xi() const
  {
    if (data["Dynamics"]["Xi"])
      return data["Dynamics"]["Xi"].as<double>();
    else if (data["Dynamics"]["BulkDriving"])
      return data["Dynamics"]["BulkDriving"].as<double>();
    else
      return defaultXi;
  }

  int Input::seed() const
  {
    if (data["Dynamics"]["Seed"])
      return data["Dynamics"]["Seed"].as<int>();
    else
      return defaultSeed;
  }
  
  double Input::nonEqAmplitude() const
  {
    if (data["Dynamics"]["NonEqForce"]) throw runtime_error("Change NonEqForce to NonEqAmplitude !");
    if (data["Dynamics"]["NonEqAmplitude"])
      return data["Dynamics"]["NonEqAmplitude"].as<double>();
    //else if (sameLetters(systemName(), "Bicolor"))
    //  return eta();
    
    return eta();
  }

  double Input::eta() const
  {
    if (data["Dynamics"]["Eta"])
      return data["Dynamics"]["Eta"].as<double>();
    else if (sameLetters(dynamicsName(), "BoundaryLangevin"))
    {
      if (data["Dynamics"]["DeltaTemperature"])
        return 2*beta() * data["Dynamics"]["DeltaTemperature"].as<double>();
      else
        return 2*beta() * defaultDeltaTemperature;
    }
    else
      return defaultEta;
  }
  
  double Input::drift() const
  {
    if (data["Dynamics"]["Drift"])
      return data["Dynamics"]["Drift"].as<double>();
    else return defaultDrift;
  }
  
  double Input::bulkDriving() const
  {
    return nu() / (2 * (nbOfParticles()-2));
    /*if (data["Dynamics"]["BulkDriving"])
      return data["Dynamics"]["BulkDriving"].as<double>();
    else 
      return nu() / (2 * (nbOfParticles()-2));*/
    //else return defaultBulkDriving;
  }
  
  double Input::nu() const
  {
    if (data["Dynamics"]["Nu"])
      return data["Dynamics"]["Nu"].as<double>();
    else if (sameLetters(dynamicsName(), "BoundaryLangevin"))
    {
      if (data["Dynamics"]["BulkDriving"])
        return (2 * (nbOfParticles()-2)) * data["Dynamics"]["BulkDriving"].as<double>();
      else
        return (2 * (nbOfParticles()-2)) * defaultBulkDriving;
    }
    else 
      return defaultNu;
  }
  
  double Input::flux() const
  {
    if (data["Dynamics"]["Flux"])
      return data["Dynamics"]["Flux"].as<double>();
    else return defaultFlux;
  }

  double Input::kappa() const
  {
    if (data["Dynamics"]["Kappa"])
      return data["Dynamics"]["Kappa"].as<double>();
    else return defaultKappa;
  }

 int Input::MTSfrequency() const
  {
    if (data["Dynamics"]["MTS frequency"])
      return data["Dynamics"]["MTS frequency"].as<int>();
    else return defaultMTSfrequency;
  }

  double Input::initialInternalTemperature() const
  {
    if (data["Dynamics"]["InitialInternalTemperature"])
      return data["Dynamics"]["InitialInternalTemperature"].as<double>();
    else return temperature();
  }

}
