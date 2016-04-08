#include "Input.hpp"

using std::cout;
using std::endl;
using std::string;

namespace simol {
  
  string Input::controlVariateName() const
  {
    if (data["ControlVariate"]["Name"])
      return data["ControlVariate"]["Name"].as<string>();
    else
      return "None";
  }

  string Input::controlVariateCoeffsPath() const
  {
    if (data["ControlVariate"]["CoeffsPath"])
      return data["ControlVariate"]["CoeffsPath"].as<string>();
    else
      assert(false);
    return "None";
  }
  
  bool Input::doGalerkinCV() const
  {
    if (data["ControlVariate"])
      if (data["Galerkin"]["Basis"] && !data["ControlVariate"]["CoeffsPath"])
	return true;
		return false;
  }
  
  bool Input::isGalerkin() const
  {
    if (data["Galerkin"])
      if (data["Galerkin"]["Resolve"])
	if (data["Galerkin"]["Resolve"].as<string>() == "Yes")
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
    
    throw std::invalid_argument("Number of Fourier modes missing");
  }
  
  int Input::nbOfHermite() const
  {
    if(data["Galerkin"])
      if(data["Galerkin"]["Basis"])
	if(data["Galerkin"]["Basis"]["Hermite"])
	  return data["Galerkin"]["Basis"]["Hermite"].as<int>();
    
    throw std::invalid_argument("Number of Hermite modes missing");
  }
  
}
