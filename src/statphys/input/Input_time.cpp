#include "simol/statphys/input/Input.hpp"

using std::cout;
using std::endl;
using std::string;

const int defaultThermalizationTime = 0;
const int defaultBurnInTime = 0;


namespace simol
{

  // timestep for the discretization
  double Input::timeStep() const
  {
    if (data["Time"]["Step"])
      return data["Time"]["Step"].as<double>();
    throw std::invalid_argument("Timestep missing");
  }

  /// Number of steps in the "output" part of the simulation
  /// Truncated to a multiple of the timeStep
  long int Input::nbOfSteps() const
  {
    static bool warningGiven = false;
    if (data["Time"]["TotalNbOfSteps"])
      return data["Time"]["TotalNbOfSteps"].as<long int>();
    else if (data["Time"]["TotalTime"])
    {
      double totalTime = data["Time"]["TotalTime"].as<double>();
      long int totalNb = totalTime / timeStep();
      double totalTrunc = totalNb * timeStep();
      if (!warningGiven && (totalTime != totalTrunc))
      {
        warningGiven = true;
        cout << "#### /!\\ totalTime truncated to " << totalTrunc << " ####" << endl;
      }
      return totalNb;
    }
    else throw std::runtime_error("Number of Iterations not specified !");
  }

  /// Number of steps in the thermalization part (temperatures imposed everywhere)
  /// Truncated to a multiple of the timeStep
  long int Input::thermalizationNbOfSteps() const
  {
    static bool warningGiven = false;
    if (data["Time"]["ThermalizationNbOfSteps"])
      return data["Time"]["ThermalizationNbOfSteps"].as<long int>();
    else if (data["Time"]["ThermalizationTime"])
    {
      double thermalTime = data["Time"]["ThermalizationTime"].as<double>();
      long int thermalNb = thermalTime / timeStep();
      double thermalTrunc = thermalNb * timeStep();
      if (!warningGiven && (thermalTime != thermalTrunc))
      {
        warningGiven = true;
        cout << "#### /!\\ thermalTime truncated to " << thermalTrunc << " ####" << endl;
      }
      return thermalNb;
    }
    else
      return defaultThermalizationTime / timeStep();
  }

  /// Number of steps in the burnIn part (no output)
  /// Truncated to a multiple of the timeStep
  long int Input::burninNbOfSteps() const
  {
    static bool warningGiven = false;
    if (data["Time"]["BurnInNbOfSteps"])
      return data["Time"]["BurnInNbOfSteps"].as<long int>();
    else if (data["Time"]["BurnInTime"])
    {
      double burninTime = data["Time"]["BurnInTime"].as<double>();
      long int burninNb = burninTime / timeStep();
      double burninTrunc = burninNb * timeStep();
      if (!warningGiven && (burninTime != burninTrunc))
      {
        warningGiven = true;
        cout << "#### /!\\ burninTime truncated to " << burninTrunc << " ####" << endl;
      }
      return burninNb;
    }
    else
      return defaultBurnInTime / timeStep();
  }

}
