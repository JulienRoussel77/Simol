#include "simol/statphys/input/Input.hpp"

using std::cout;
using std::endl;
using std::string;

namespace simol
{
  string Input::CVObservable() const
  {
    if (data["ControlVariate"])
      if (data["ControlVariate"]["Observable"])
        return data["ControlVariate"]["Observable"].as<string>();

    return "None";
  }
  
  bool Input::doCVObservable(int idObs) const
  {
    //cout << "doCVObservable : " << nameOfObs(idObs) << " =?= " << CVObservable() << endl;
    return nameOfObs(idObs) == CVObservable();
  }

  string Input::controlVariateName() const
  {
    if (data["ControlVariate"])
      if (data["ControlVariate"]["Name"])
        return data["ControlVariate"]["Name"].as<string>();

    return "None";
  }
  
  string Input::galerkinPotentialName() const
  {
    if (data["Galerkin"])
      if (data["Galerkin"]["PotentialName"])
        return data["Galerkin"]["PotentialName"].as<string>();

    if (externalPotentialName() == "Harmonic")
      return "HarmonicFE";
    if (externalPotentialName() == "DoubleWell")
      return "DoubleWellFE";
    
    return externalPotentialName();
  }

  string Input::controlVariateCoeffsPath() const
  {
    if (data["ControlVariate"])
      if (data["ControlVariate"]["CoeffsPath"])
        return data["ControlVariate"]["CoeffsPath"].as<string>();

    return "None";
  }

  bool Input::doGalerkinCV() const
  {
    //if (data["ControlVariate"])
      //if (data["Galerkin"]["Basis"] && !data["ControlVariate"]["CoeffsPath"])
    if (data["Galerkin"])
        return true;
    return false;
  }

  bool Input::isGalerkin() const
  {
    if (data["Galerkin"])
      if (data["Galerkin"]["Resolve"])
        if (isYes(data["Galerkin"]["Resolve"].as<string>()))
          return true;

    return false;
  }

  string Input::basisElts() const
  {
    if (data["Basis"])
      if (data["Basis"]["Elements"])
        return data["Basis"]["Elements"].as<string>();
    return "None";
  }

  int Input::nbOfQModes() const
  {
    if (data["Basis"])
        if (data["Basis"]["QModes"])
          return data["Basis"]["QModes"].as<int>();

    return 0;
    //throw std::invalid_argument("Number of Fourier modes missing");
  }

  int Input::nbOfPModes() const
  {
    if(data["Basis"])
      if(data["Basis"]["PModes"])
        return data["Basis"]["PModes"].as<int>();

    return 0;
    //throw std::invalid_argument("Number of Hermite modes missing");
  }
  
  bool Input::doGalerkinNonequilibrium() const
  {
    if(data["Galerkin"])
      if(data["Galerkin"]["Nonequilibrium"])
        if (isYes(data["Galerkin"]["Nonequilibrium"].as<string>()))
          return true;
    return false;
    //throw std::invalid_argument("Number of Hermite modes missing");
  }
  
  bool Input::doComputeRef() const
  {
    if(data["Galerkin"])
      if(data["Galerkin"]["ComputeRef"])
        if (isYes(data["Galerkin"]["ComputeRef"].as<string>()))
          return true;
    return false;
    //throw std::invalid_argument("Number of Hermite modes missing");
  }
  
  bool Input::fitModiFlow() const
  {
    if (data["ControlVariate"])
      if (data["ControlVariate"]["fitModiFlow"])
        if (isYes(data["Galerkin"]["fitModiFlow"].as<string>()))
          return true;
    return false;
  }
  
  string Input::CVOperator() const
  {
    if (data["ControlVariate"])
      if (data["ControlVariate"]["Operator"])
        return data["ControlVariate"]["Operator"].as<string>();

    return dynamicsName();
  }
  
  double Input::omegaHermite() const
  {
    if (data["Basis"])
      if (data["Basis"]["Omega"])
        return data["Basis"]["Omega"].as<double>();

    return 1;
    //throw std::invalid_argument("Number of Fourier modes missing");
  }
  
  double Input::meshStep() const
  {
    if (data["Basis"])
      if (data["Basis"]["MeshStep"])
        return data["Basis"]["MeshStep"].as<double>();

    return 0;
    //throw std::invalid_argument("Number of Fourier modes missing");
  }
  
  double Input::integrationStep() const
  {
    if (data["Galerkin"])
      if (data["Galerkin"]["Integration"])
        if (data["Galerkin"]["Integration"]["Step"])
          return data["Galerkin"]["Integration"]["Step"].as<double>();

    return 0.01;
  }
    
  double Input::integrationQMin() const
  {
    if (data["Galerkin"])
      if (data["Galerkin"]["Integration"])
        if (data["Galerkin"]["Integration"]["QMin"])
          return data["Galerkin"]["Integration"]["QMin"].as<double>();

    return 0;
  }
  
  double Input::integrationLength() const
  {
    if (data["Galerkin"])
      if (data["Galerkin"]["Integration"])
        if (data["Galerkin"]["Integration"]["Length"])
          return data["Galerkin"]["Integration"]["Length"].as<double>();

    return 0;
  }
  
  double Input::galerkinEpsilon() const
  {
    if (data["Galerkin"])
      if (data["Galerkin"]["Epsilon"])
          return data["Galerkin"]["Epsilon"].as<double>();

    return 0;
  }
}
