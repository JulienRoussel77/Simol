#include "simol/statphys/input/Input.hpp"

using std::cout;
using std::endl;
using std::string;

const double defaultSeed = 0;
const double defaultGamma = 1;
const double defaultXi = 0;
const double defaultDeltaTemperature = 0;
const double defaultExternalForce = 0;
const double defaultTauBending = 0;
const double defaultHeatCapacity = 1;

namespace simol
{

  string Input::dynamicsName() const {return data["Dynamics"]["Name"].as<string>();}

  double Input::gamma() const
  {
    if (data["Dynamics"]["Gamma"])
      return data["Dynamics"]["Gamma"].as<double>();
    else throw std::runtime_error("No Gamma in the input file !");
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
    else throw std::runtime_error("Temperature not precised !");
  }

  /*double Input::temperatureLeft() const
  {
    if (data["Dynamics"]["TemperatureLeft"])
      return data["Dynamics"]["TemperatureLeft"].as<double>();
    else throw std::runtime_error("No TemperatureLeft in the input file !");
  }

  double Input::temperatureRight() const
  {
    if (data["Dynamics"]["TemperatureRight"])
      return data["Dynamics"]["TemperatureRight"].as<double>();
    else throw std::runtime_error("No TemperatureRight in the input file !");
  }*/

  double Input::beta() const
  {
    return 1 / temperature();
  }

  /*double Input::betaLeft() const
  {
    if (data["Dynamics"]["BetaLeft"])
      return data["Dynamics"]["BetaLeft"].as<double>();
    else if (data["Dynamics"]["TemperatureLeft"])
      return 1 / data["Dynamics"]["TemperatureLeft"].as<double>();
    else throw std::runtime_error("No BetaLeft in the input file !");
  }

  double Input::betaRight() const
  {
    if (data["Dynamics"]["BetaRight"])
      return data["Dynamics"]["BetaRight"].as<double>();
    else if (data["Dynamics"]["TemperatureRight"])
      return 1 / data["Dynamics"]["TemperatureRight"].as<double>();
    else throw std::runtime_error("No BetaRight in the input file !");
  }*/

  double Input::deltaTemperature() const
  {
    if (data["Dynamics"]["DeltaTemperature"])
      return data["Dynamics"]["DeltaTemperature"].as<double>();
    else return defaultDeltaTemperature;
  }

  double Input::externalForce() const
  {
    if (data["Dynamics"]["Force"])
      return data["Dynamics"]["Force"].as<double>();
    else return defaultExternalForce;
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

  double Input::eta() const
  {
    if (dynamicsName() == "BoundaryLangevin")
      return deltaTemperature();
    else
      return externalForce();
  }

  double Input::heatCapacity() const
  {
    if (data["Physics"]["Model"]["HeatCapacity"])
      return data["Physics"]["Model"]["HeatCapacity"].as<double>();
    else
      return defaultHeatCapacity;
  }


}
