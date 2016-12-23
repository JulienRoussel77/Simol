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
    return nameOfObs(idObs) == CVObservable();
  }

  string Input::controlVariateName() const
  {
    if (data["ControlVariate"])
      if (data["ControlVariate"]["Name"])
        return data["ControlVariate"]["Name"].as<string>();

    return "None";
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
        if (sameLetters(data["Galerkin"]["Resolve"].as<string>(), "yes"))
          return true;

    return false;
  }

  string Input::galerkinElts() const
  {
    if (data["Galerkin"])
      if (data["Galerkin"]["Basis"])
        if (data["Galerkin"]["Basis"]["Elements"])
          return data["Galerkin"]["Basis"]["Elements"].as<string>();
    return "None";
  }

  int Input::nbOfFourier() const
  {
    if (data["Galerkin"])
      if (data["Galerkin"]["Basis"])
        if (data["Galerkin"]["Basis"]["Fourier"])
          return data["Galerkin"]["Basis"]["Fourier"].as<int>();

    return 0;
    //throw std::invalid_argument("Number of Fourier modes missing");
  }

  int Input::nbOfHermite() const
  {
    if(data["Galerkin"])
      if(data["Galerkin"]["Basis"])
        if(data["Galerkin"]["Basis"]["Hermite"])
          return data["Galerkin"]["Basis"]["Hermite"].as<int>();

    return 0;
    //throw std::invalid_argument("Number of Hermite modes missing");
  }
  
  bool Input::doGalerkinNonequilibrium() const
  {
    if(data["Galerkin"])
      if(data["Galerkin"]["Nonequilibrium"])
        if (sameLetters(data["Galerkin"]["Nonequilibrium"].as<string>(), "yes"))
          return true;
    return false;
    //throw std::invalid_argument("Number of Hermite modes missing");
  }
  
  bool Input::doComputeRef() const
  {
    if(data["Galerkin"])
      if(data["Galerkin"]["ComputeRef"])
        if (sameLetters(data["Galerkin"]["ComputeRef"].as<string>(), "yes"))
          return true;
    return false;
    //throw std::invalid_argument("Number of Hermite modes missing");
  }
  
}
