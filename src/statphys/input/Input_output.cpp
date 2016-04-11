#include "Input.hpp"

using std::cout;
using std::endl;
using std::string;
using std::min;
using std::max;

const int defaultOutputPeriod = 1;         
const int defaultOutputProfilePeriod = 1;  
const int maxNbOfAutocoPts = 1000;

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
      name += data["Output"]["FolderName"].as<string>()+"/";
    
    return name;
  }
  

  int Input::decorrelationNbOfIterations() const
  {
    if (data["Output"]["DecorrelationTime"])
      return data["Output"]["DecorrelationTime"].as<double>() / timeStep();
    else
      return decorrelationTime() / timeStep();
  }

  double Input::decorrelationTime() const
  {
    if (data["Output"]["DecorrelationTime"])
      return data["Output"]["DecorrelationTime"].as<double>();
    else return 3*nbOfParticles();
  }

  int Input::outputPeriodNbOfIterations() const {
    if (data["Output"]["Period"])
      return data["Output"]["Period"].as<double>() / timeStep();
    else
      return defaultOutputPeriod;
  }

  double Input::outputPeriodTime() const {
    if (data["Output"]["Period"])
      return data["Output"]["Period"].as<double>();
    else
      return outputPeriodNbOfIterations() * timeStep();
  }

  int Input::outputProfilePeriodNbOfIterations() const {
    if (data["Output"]["ProfilePeriod"])
      return data["Output"]["ProfilePeriod"].as<double>() / timeStep();
    else
      return defaultOutputProfilePeriod;
  }

  double Input::outputProfilePeriodTime() const {
    if (data["Output"]["ProfilePeriod"])
      return data["Output"]["ProfilePeriod"].as<double>();
    else
      return outputProfilePeriodNbOfIterations() * timeStep();
  }

  /// Contains the number of values in an autocorrelation ProfilePeriod
  /// /!\ Causes a memory crash for the larger chains if too big
  int Input::nbOfAutocoPts() const
  {
    return min(maxNbOfAutocoPts, (int)decorrelationNbOfIterations());
  }
  
  bool Input::doFinalFlow() const
  {
    if (data["Output"]["doFinalFlow"])
      if (data["Output"]["doFinalFlow"].as<string>() == "no")
	return false;
    return true;
  }
  
  bool Input::doFinalVelocity() const
  {
    if (data["Output"]["doFinalVelocity"])
      if (data["Output"]["doFinalVelocity"].as<string>() == "no")
	return false;
    return true;
  }
  
}
