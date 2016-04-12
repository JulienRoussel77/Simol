#include "Input.hpp"

using std::cout;
using std::endl;
using std::string;
using std::min;
using std::max;

const int defaultPrintPeriod = 1;         
const int defaultLongPrintPeriod = 1;  
const int maxNbOfAutocoPts = 1000;
const int defaultDecorrelationTime = 0; 

namespace simol {

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
      return "../../../output/"+dynamicsName()+"/"+systemName()+"/"+potentialName()+"/";
    else 
      return "../../../output/";
  }
  
  string Input::parametersName() const
  {
    string name = simuTypeName();
        
    if (controlVariateName() != "None")
      name += controlVariateName()+"/";
    
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
  

  int Input::decorrelationNbOfSteps() const
  {
    if (data["Output"]["DecorrelationTime"])
      return data["Output"]["DecorrelationTime"].as<double>() / timeStep();
    else
      return defaultDecorrelationTime;
  }

  double Input::decorrelationTime() const
  {
    if (data["Output"]["DecorrelationTime"])
      return data["Output"]["DecorrelationTime"].as<double>();
    else 
      return defaultDecorrelationTime;
  }

  int Input::printPeriodNbOfSteps() const {
    if (data["Output"]["PrintPeriodNbOfSteps"])
      return data["Output"]["PrintPeriodNbOfSteps"].as<int>();
    else if (data["Output"]["PrintPeriod"])
    {
      if (modulo(data["Output"]["PrintPeriod"].as<double>(), 0, timeStep()) != 0)
        cout << "#### PrintPeriod truncated ####" << endl;
      return data["Output"]["PrintPeriod"].as<double>() / timeStep();
    }
    else
      return defaultPrintPeriod;
  }

  double Input::printPeriodTime() const {
      return printPeriodNbOfSteps() * timeStep();
  }

  int Input::printLongPeriodNbOfSteps() const {
    if (data["Output"]["LongPrintPeriodNbOfSteps"])
      return data["Output"]["LongPrintPeriodNbOfSteps"].as<int>();
    if (data["Output"]["LongPrintPeriod"])
    {
      if (modulo(data["Output"]["LongPrintPeriod"].as<double>(), 0, timeStep()) != 0)
        cout << "#### LongPrintPeriod truncated ####" << endl;
        return data["Output"]["LongPrintPeriod"].as<double>() / timeStep();
    }
    else
      return defaultLongPrintPeriod;
  }

  double Input::printLongPeriodTime() const {
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
