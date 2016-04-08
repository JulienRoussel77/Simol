#include "Input.hpp"

using std::cout;
using std::endl;
using std::string;

namespace simol {

  // timestep for the discretization 
  double Input::timeStep() const
  {
    if (data["Time"]["Step"])
      return data["Time"]["Step"].as<double>();
    throw std::invalid_argument("Timestep missing");
  }
  
  // Number of iterations in the "output" part of the simulation
  int Input::nbOfIterations() const
  {
    if (data["Time"]["Number"])
      return data["Time"]["Number"].as<int>();
    else if (data["Time"]["FinalTime"])
      return data["Time"]["FinalTime"].as<double>() / timeStep();
    else
      {cout << "Nb of Iterations not specified !" << endl;exit(1);}
  }
  
  // Number of iterations in the thermalization part (temperatures imposed everywhere)
  int Input::nbOfThermalIterations() const
  {
    if (data["Time"]["ThermalNumber"])
      return data["Time"]["ThermalNumber"].as<int>();
    else if (data["Time"]["ThermalTime"])
      return data["Time"]["ThermalTime"].as<double>() / timeStep();
    else
      return 0;
  }

  // Number of iterations in the burning part (no output)
  int Input::nbOfBurningIterations() const
  {
    if (data["Time"]["BurningNumber"])
      return data["Time"]["BurningNumber"].as<int>();
    else if (data["Time"]["BurningTime"])
      return data["Time"]["BurningTime"].as<double>() / timeStep();
    else
      return 0;
  }

}
