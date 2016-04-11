#include "Input.hpp"

using std::cout;
using std::endl;
using std::string;

const int defaultThermalIterations = 0;
const int defaultBurnInIterations = 0;


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
    else throw std::runtime_error("Nb of Iterations not specified !");
  }
  
  // Number of iterations in the thermalization part (temperatures imposed everywhere)
  int Input::nbOfThermalIterations() const
  {
    if (data["Time"]["ThermalNumber"])
      return data["Time"]["ThermalNumber"].as<int>();
    else if (data["Time"]["ThermalTime"])
      return data["Time"]["ThermalTime"].as<double>() / timeStep();
    else
      return defaultThermalIterations;
  }

  // Number of iterations in the burnIn part (no output)
  int Input::nbOfBurnInIterations() const
  {
    if (data["Time"]["BurnInNumber"])
      return data["Time"]["BurnInNumber"].as<int>();
    else if (data["Time"]["BurnInTime"])
      return data["Time"]["BurnInTime"].as<double>() / timeStep();
    else
      return defaultBurnInIterations;
  }

}
