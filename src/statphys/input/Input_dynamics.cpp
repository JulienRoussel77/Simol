#include "Input.hpp"

using std::cout;
using std::endl;
using std::string;

const double defaultSeed = 0;
const double defaultGamma = 1; 
const double defaultXi = 0;
const double defaultExternalForce = 0;
const double defaultTauBending = 0;
const double defaultHeatCapacity = 1;

namespace simol {

  string Input::dynamicsName() const {return data["Dynamics"]["Name"].as<string>();}
  
  double Input::gamma() const {return data["Dynamics"]["Gamma"].as<double>();}
  
  double Input::temperature() const
  {
    if (data["Dynamics"]["Temperature"])
      return data["Dynamics"]["Temperature"].as<double>();
    else if (data["Dynamics"]["Beta"])
      return 1 / data["Dynamics"]["Beta"].as<double>();
    else if (data["Dynamics"]["TemperatureLeft"] && data["Dynamics"]["TemperatureRight"])
      return (temperatureLeft() + temperatureRight()) / 2;
    else if (data["Dynamics"]["BetaLeft"] && data["Dynamics"]["BetaRight"])
      return .5/betaLeft() + .5/betaRight();
    else
      {cout << "Temperature not precised !" << endl;exit(1);}
  }

  double Input::temperatureLeft() const
  {
    return data["Dynamics"]["TemperatureLeft"].as<double>();
  }

  double Input::temperatureRight() const
  {
    return data["Dynamics"]["TemperatureRight"].as<double>();
  }

  double Input::beta() const
  {
    if (data["Dynamics"]["Beta"])
      return data["Dynamics"]["Beta"].as<double>();
    else if (data["Dynamics"]["Temperature"])
      return 1 / data["Dynamics"]["Temperature"].as<double>();
    else if (data["Dynamics"]["BetaLeft"] && data["Dynamics"]["BetaRight"])
      return 2. / (1/betaLeft() + 1/betaRight());
    else if (data["Dynamics"]["TemperatureLeft"] && data["Dynamics"]["TemperatureRight"])
      return 2/(temperatureLeft() + temperatureRight());
    else
      throw std::invalid_argument("Beta not precised !");
  }

  double Input::betaLeft() const
  {
    if (data["Dynamics"]["BetaLeft"])
      return data["Dynamics"]["BetaLeft"].as<double>();
    else if (data["Dynamics"]["TemperatureLeft"])
      return 1 / data["Dynamics"]["TemperatureLeft"].as<double>();
    else assert(false);
  }

  double Input::betaRight() const
  {
    if (data["Dynamics"]["BetaRight"])
      return data["Dynamics"]["BetaRight"].as<double>();
    else if (data["Dynamics"]["TemperatureRight"])
      return 1 / data["Dynamics"]["TemperatureRight"].as<double>();
    else assert(false);
  }

  double Input::externalForce() const {
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
      return (temperatureLeft() - temperatureRight())/2;
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
