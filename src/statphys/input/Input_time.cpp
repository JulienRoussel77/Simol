#include "Input.hpp"

using std::cout;
using std::endl;
using std::string;

const int defaultThermalizationNbOfSteps = 0;
const int defaultBurnInNbOfSteps = 0;


namespace simol {

  // timestep for the discretization 
  double Input::timeStep() const
  {
    if (data["Time"]["Step"])
      return data["Time"]["Step"].as<double>();
    throw std::invalid_argument("Timestep missing");
  }
  
  // Number of steps in the "output" part of the simulation
  int Input::nbOfSteps() const
  {
    if (data["Time"]["TotalNbOfSteps"])
      return data["Time"]["TotalNbOfSteps"].as<int>();
    else if (data["Time"]["TotalTime"])
      return data["Time"]["TotalTime"].as<double>() / timeStep();
    else throw std::runtime_error("Number of Iterations not specified !");
  }
  
  // Number of steps in the thermalization part (temperatures imposed everywhere)
  int Input::thermalizationNbOfSteps() const
  {
    if (data["Time"]["ThermalizationNbOfSteps"])
      return data["Time"]["ThermalizationNbOfSteps"].as<int>();
    else if (data["Time"]["ThermalizationTime"])
      return data["Time"]["ThermalizationTime"].as<double>() / timeStep();
    else
      return defaultThermalizationNbOfSteps;
  }

  // Number of steps in the burnIn part (no output)
  int Input::burninNbOfSteps() const
  {
    if (data["Time"]["BurnInNbOfSteps"])
      return data["Time"]["BurnInNbOfSteps"].as<int>();
    else if (data["Time"]["BurnInTime"])
      return data["Time"]["BurnInTime"].as<double>() / timeStep();
    else
      return defaultBurnInNbOfSteps;
  }

}
