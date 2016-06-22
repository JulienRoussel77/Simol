#include "simol/statphys/input/Input.hpp"

using std::cout;
using std::endl;
using std::string;
using std::min;
using std::max;

const int defaultPrintPeriod = 0;
const int defaultLongPrintPeriod = 0;
const int maxNbOfAutocoPts = 1000;
const int defaultDecorrelationTime = 0;

namespace simol
{

  ///Transforms a double into a nice string
  /// 1.30000 -> 1.3     1.000 -> 1
  string doubleToString(double x)
  {
    string s = to_string (x);
    if (s.find(".") != std::string::npos)
      s.erase ( s.find_last_not_of('0') + 1, string::npos );
    s.erase ( s.find_last_not_of('.') + 1, string::npos );
    return s;
  }

  string Input::simuTypeName() const
  {
    if (data["Output"]["SimuTypeName"]
        && data["Output"]["SimuTypeName"].as<string>() == "yes")
      //return "../../../output/" + dynamicsName() + "/" + systemName() + "/" + potentialName() + "/";
      return "output/" + dynamicsName() + "/" + systemName() + "/" + potentialName() + "/";
    else
      //return "../../../output/";
      return "output/";
  }

  string Input::parametersName() const
  {
    string name = simuTypeName();

    if (controlVariateName() != "None")
      name += controlVariateName() + "/";

    if (data["Output"]["ParametersName"]
        && data["Output"]["ParametersName"].as<string>() == "yes")
    {
      if (dynamicsName() == "BoundaryLangevin")
        name += "N" + to_string(nbOfParticles()) + "_";

      name += "dt" + doubleToString(timeStep()) + "_eta" + doubleToString(eta());

      if (dynamicsName() == "BoundaryLangevin")
        name += "_xi" + doubleToString(xi());

      name += "/";
    }
    return name;
  }

  string Input::outputFolderName() const
  {
    string name = parametersName();

    if (data["Output"]["FolderName"])
      name += data["Output"]["FolderName"].as<string>() + "/";

    return name;
  }

  ///
  ///Returns the a number of steps
  int Input::decorrelationNbOfSteps() const
  {
    static bool warningGiven = false;
    if (data["Output"]["DecorrelationTime"])
    {
      double decoTime = data["Output"]["DecorrelationTime"].as<double>();  // decorrelation time not truncated
      int decoNb = (int)(decoTime / timeStep());
      double decoTrunc = decoNb * timeStep();                 // decorrelation time truncated
      if (!warningGiven && (decoTime != decoTrunc))           //sends a warning the first time the function is called
      {
        warningGiven = true;
        cout << "#### /!\\ DecorrelationTime truncated to " << decoTrunc << " ####" << endl;
      }
      return decoNb;
    }
    else
      return (int)(defaultDecorrelationTime / timeStep());

  }
  ///
  ///Returns the truncated decorrelation time, which is a multiple of the one given in input
  double Input::decorrelationTime() const
  {
    return decorrelationNbOfSteps() * timeStep();
  }

  int Input::printPeriodNbOfSteps() const
  {
    static bool warningGiven = false;
    if (data["Output"]["PrintPeriodNbOfSteps"])
      return data["Output"]["PrintPeriodNbOfSteps"].as<int>();
    else if (data["Output"]["PrintPeriod"])
    {
      double periodTime = data["Output"]["PrintPeriod"].as<double>();
      int periodNb = (int)( periodTime / timeStep());
      double periodTrunc = periodNb * timeStep();
      if (!warningGiven && (periodTime != periodTrunc))
      {
        warningGiven = true;
        cout << "#### /!\\ PrintPeriod truncated to " << periodTrunc << " ####" << endl;
      }
      return periodNb;
    }
    else
      return defaultPrintPeriod / timeStep();
  }

  double Input::printPeriodTime() const
  {
    return printPeriodNbOfSteps() * timeStep();
  }

  int Input::printLongPeriodNbOfSteps() const
  {
    static bool warningGiven = false;
    if (data["Output"]["LongPrintPeriodNbOfSteps"])
      return data["Output"]["LongPrintPeriodNbOfSteps"].as<int>();
    if (data["Output"]["LongPrintPeriod"])
    {
      double periodTime = data["Output"]["LongPrintPeriod"].as<double>();
      int periodNb = (int)( periodTime / timeStep());
      double periodTrunc = periodNb * timeStep();
      if (!warningGiven && (periodTime != periodTrunc))
      {
        warningGiven = true;
        cout << "#### /!\\ LongPrintPeriod truncated to " << periodTrunc << " ####" << endl;
      }
      return periodTrunc;
    }
    else
      return defaultLongPrintPeriod / timeStep();
  }

  double Input::printLongPeriodTime() const
  {
    return printLongPeriodNbOfSteps() * timeStep();
  }

  /// Contains the number of values in an autocorrelation LongPeriod
  /// /!\ Causes a memory crash for the larger chains if too big
  int Input::nbOfAutocoPts() const
  {
    return min(maxNbOfAutocoPts, (int)decorrelationNbOfSteps());
  }

  bool Input::doFinalFlow() const
  {
    if (data["Output"]["doFinalFlow"])
      if (data["Output"]["doFinalFlow"].as<string>() == "yes")
        return true;
    return false;
  }

  bool Input::doFinalVelocity() const
  {
    if (data["Output"]["doFinalVelocity"])
      if (data["Output"]["doFinalVelocity"].as<string>() == "yes")
        return true;
    return false;
  }

}
