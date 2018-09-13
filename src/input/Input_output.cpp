#include "simol/input/Input.hpp"

using std::cout;
using std::endl;
using std::string;
using std::min;
using std::max;

const int defaultPrintPeriod = 0;
const int defaultPrintLongPeriod = 0;
const int maxNbOfAutocoPts = 200;
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
    string name = "output/";
    if (data["Output"]["SimuTypeName"]
        && sameLetters(data["Output"]["SimuTypeName"].as<string>(), "yes"))
    {
      name += dynamicsName() + "/" + systemName() + "/";
      if (potentialName() == "TwoTypes")
        name += secondPotentialName() + "/" + firstPotentialName() + "/";
      else
      name += potentialName() + "/";
    }
    
    if (data["Output"]["FolderName"])
      name += data["Output"]["FolderName"].as<string>() + "/";

    //if (controlVariateName() != "None")
    //  name += controlVariateName() + "/";
    
    return name;
    
    
    /*if (data["Output"]["SimuTypeName"]
        && sameLetters(data["Output"]["SimuTypeName"].as<string>(), "yes"))
      //return "../../../output/" + dynamicsName() + "/" + systemName() + "/" + potentialName() + "/";
      return "output/" + dynamicsName() + "/" + systemName() + "/" + potentialName() + "/";
    else
      //return "../../../output/";
      return "output/";*/
  }

  string Input::outputFolderName() const
  {
    string name = simuTypeName();
    
    if (data["Output"]["ParametersName"])
    {
      //YAML does not seem to detect automatically the type (scalar or array) of the input, so we have to test it manually
      
      vector<string> vecParameters;
      if (data["Output"]["ParametersName"].Type() == YAML::NodeType::Scalar)
        vecParameters = vector<string>(1, data["Output"]["ParametersName"].as<string>());
      else if (data["Output"]["ParametersName"].Type() == YAML::NodeType::Sequence)
      //for (int iOfVec = 0; iOfVec < (int)data["Output"]["ParametersName"].size(); iOfVec++)
      vecParameters = data["Output"]["ParametersName"].as<std::vector<string>>();
      
      for (int iOfVec = 0; iOfVec < (int)vecParameters.size(); iOfVec++)
      {
        if (iOfVec != 0) name += "_";
        name += vecParameters[iOfVec];
        if (vecParameters[iOfVec] == "dt")            name += doubleToString(timeStep());
        else if (vecParameters[iOfVec] == "N")        name += doubleToString(nbOfParticles());
        else if (vecParameters[iOfVec] == "eta")      name += doubleToString(eta());
        else if (vecParameters[iOfVec] == "xi")       name += doubleToString(xi());
        else if (vecParameters[iOfVec] == "beta")     name += doubleToString(beta());
        else if (vecParameters[iOfVec] == "T")        name += doubleToString(temperature());
        else if (vecParameters[iOfVec] == "gamma")    name += doubleToString(gamma());
        else if (vecParameters[iOfVec] == "seed")     name += doubleToString(seed());
        else if (vecParameters[iOfVec] == "potAlpha") name += doubleToString(potentialAlpha());
        else if (vecParameters[iOfVec] == "potBeta")  name += doubleToString(potentialBeta());
        else if (vecParameters[iOfVec] == "c")  name += doubleToString(potentialCenter());
        else if (vecParameters[iOfVec] == "MTSfrequency")  name += doubleToString(MTSfrequency());
        else if (vecParameters[iOfVec] == "drift")  name += doubleToString(drift());
        else if (vecParameters[iOfVec] == "a")  name += doubleToString(latticeParameter());
        else if (vecParameters[iOfVec] == "sol")  name += doubleToString(interactionRatio());
        else if (vecParameters[iOfVec] == "bulk")  name += doubleToString(bulkDriving());
        else if (vecParameters[iOfVec] == "DT")  name += doubleToString(deltaTemperature());
        else if (vecParameters[iOfVec] == "flux")  name += doubleToString(flux());
        else if (vecParameters[iOfVec] == "delta")  name += doubleToString(coupling());
        else if (vecParameters[iOfVec] == "nu")  name += doubleToString(nu());
        else throw std::runtime_error(vecParameters[iOfVec] + " is not a parameter name !");
      }
      name += "/";
    }
    return name;
  }

  ///
  ///Returns the a number of steps
  int Input::decorrelationNbOfSteps() const
  {
    // This static bools prevents from giving several times the same warning
    static bool warningGiven = false;
    double decoTime;
    if (data["Output"]["DecorrelationTime"]) decoTime = data["Output"]["DecorrelationTime"].as<double>();  // decorrelation time not truncated
    else decoTime = defaultDecorrelationTime;
    int decoNb = (int)(decoTime / timeStep());
    double decoTrunc = decoNb * timeStep();                 // decorrelation time truncated
    if (!warningGiven && (decoTime != decoTrunc))           //sends a warning the first time the function is called
    {
      warningGiven = true;
      cout << "#### /!\\ DecorrelationTime truncated to " << decoTrunc << " ####" << endl;
    }
    return decoNb;
  }
  ///
  ///Returns the truncated decorrelation time, which is a multiple of the one given in input
  double Input::decorrelationTime() const
  {
    return decorrelationNbOfSteps() * timeStep();
  }
  
  ///
  ///Returns the a number of steps
  int Input::shortDecorrelationNbOfSteps() const
  {
    // This static bools prevents from giving several times the same warning
    static bool warningGiven = false;
    double decoTime;
    if (data["Output"]["ShortDecorrelationTime"]) decoTime = data["Output"]["ShortDecorrelationTime"].as<double>();  // decorrelation time not truncated
    else decoTime = decorrelationTime();
    int decoNb = (int)(decoTime / timeStep());
    double decoTrunc = decoNb * timeStep();                 // decorrelation time truncated
    if (!warningGiven && (decoTime != decoTrunc))           //sends a warning the first time the function is called
    {
      warningGiven = true;
      cout << "#### /!\\ ShortDecorrelationTime truncated to " << decoTrunc << " ####" << endl;
    }
    return decoNb;
  }
  ///
  ///Returns the truncated decorrelation time, which is a multiple of the one given in input
  double Input::shortDecorrelationTime() const
  {
    return shortDecorrelationNbOfSteps() * timeStep();
  }

  int Input::printPeriodNbOfSteps() const
  {
    // This static bools prevents from giving several times the same warning
    static bool warningGiven = false;
    if (data["Output"]["PrintPeriodNbOfSteps"])
      return data["Output"]["PrintPeriodNbOfSteps"].as<int>();
    
    double periodTime;
    if (data["Output"]["PrintPeriod"]) periodTime = data["Output"]["PrintPeriod"].as<double>();
    else periodTime = defaultPrintPeriod;
    int periodNb = ceil( periodTime / timeStep());
    double periodTrunc = periodNb * timeStep();
    if (!warningGiven && (periodTime != periodTrunc))
    {
      warningGiven = true;
      cout << "#### /!\\ PrintPeriod truncated to " << periodTrunc << " ####" << endl;
    }
    return periodNb;
  }

  double Input::printPeriodTime() const
  {
    return printPeriodNbOfSteps() * timeStep();
  }

  int Input::printLongPeriodNbOfSteps() const
  {
    // This static bools prevents from giving several times the same warning
    static bool warningGiven = false;
    if (data["Output"]["PrintLongPeriodNbOfSteps"])
      return data["Output"]["PrintLongPeriodNbOfSteps"].as<int>();
    
    double periodTime;
    if (data["Output"]["PrintLongPeriod"]) periodTime = data["Output"]["PrintLongPeriod"].as<double>();
    else periodTime = defaultPrintLongPeriod;
    int periodNb = ceil( periodTime / timeStep());
    double periodTrunc = periodNb * timeStep();
    if (!warningGiven && (periodTime != periodTrunc))
    {
      warningGiven = true;
      cout << "#### /!\\ PrintLongPeriod truncated to " << periodTrunc << " ####" << endl;
    }
    return periodNb;
  }

  double Input::printLongPeriodTime() const
  {
    return printLongPeriodNbOfSteps() * timeStep();
  }
  
  ///
  /// Number of successive iterations that are gathered for an autocorrelation point
  /// If maxNbOfAutocoPts is not reached, returns 1
  int Input::autocoPtsBinSize() const
  {
    return max(1, (int)decorrelationNbOfSteps()/maxNbOfAutocoPts);
  }

  /// Contains the number of values in an autocorrelation LongPeriod
  /// /!\ Causes a memory crash for the larger chains if too big
  int Input::nbOfAutocoPts() const
  {
    return decorrelationNbOfSteps() / autocoPtsBinSize();
    //return min(maxNbOfAutocoPts, (int)decorrelationNbOfSteps());
  }
  
  ///
  /// Number of successive iterations that are gathered for an autocorrelation point
  /// If maxNbOfAutocoPts is not reached, returns 1
  int Input::shortAutocoPtsBinSize() const
  {
    return max(1, (int)shortDecorrelationNbOfSteps()/maxNbOfAutocoPts);
  }
  
  /// Contains the number of values in an autocorrelation LongPeriod
  /// /!\ Causes a memory crash for the larger chains if too big
  int Input::nbOfShortAutocoPts() const
  {
    return shortDecorrelationNbOfSteps() / shortAutocoPtsBinSize();
    //return min(maxNbOfAutocoPts, (int)shortDecorrelationNbOfSteps());
  }
  
  
  
    
  //---------------- check whether observables should be add --------------------
  
  string Input::nameOfObs(int idObs) const
  {
    switch (idObs)
    {
      case idKineticEnergy : return "kineticEnergy";
      case idPotentialEnergy : return "potentialEnergy";
      case idTotalEnergy : return "totalEnergy";
      case idPressure : return "pressure";
      case idInternalEnergy : return "internalEnergy";
      case idInternalTemperature : return "internalTemperature";
      case idVelocity : return "velocity";
      case idForce : return "force";
      case idLength : return "length";
      case idMidFlux : return "midFlux";
      case idSumFlux : return "sumFlux";
      case idModiFlux : return "modiFlux";
      case idLagrangeMultiplier: return "lagrangeMultiplier";
      default : throw std::runtime_error("This observable id corresponds to no observable !");
    }
  }
  
  bool Input::doObservable(int idObs) const
  {
    switch (idObs)
    {
      case idKineticEnergy : return doKineticEnergy();
      case idPotentialEnergy : return doPotentialEnergy();
      case idTotalEnergy : return doTotalEnergy();
      case idPressure : return doPressure();
      case idInternalEnergy : return doInternalEnergy();
      case idInternalTemperature : return doInternalTemperature();
      case idVelocity : return doVelocity();
      case idForce : return doForce();
      case idLength : return doLength();
      case idMidFlux : return doMidFlux();
      case idSumFlux : return doSumFlux();
      case idModiFlux : return doModiFlux();
      case idLagrangeMultiplier: return doLagrangeMultiplier();
      default : throw std::runtime_error("doObs not defined for this observable id !");
    }
  }
      
  bool Input::doKineticEnergy() const
  {return dynamicsName() != "Overdamped";}
  
  bool Input::doPotentialEnergy() const
  {return true;}
  
  bool Input::doTotalEnergy() const
  {return true;}
  
  bool Input::doPressure() const
  {return true;}
  
  bool Input::doDPDE() const
  {return dynamicsName() == "DPDE";}
  
  bool Input::doInternalEnergy() const
  {return dynamicsName() == "DPDE";}
  
  bool Input::doInternalTemperature() const
  {return dynamicsName() == "DPDE";}
  
  bool Input::doVelocity() const
  {return (dynamicsName() != "Overdamped" && (systemName() == "Isolated" || systemName() == "Bicolor"));}
  
  bool Input::doForce() const
  {return (dynamicsName() == "Overdamped" || systemName() == "Bicolor");}
  
  bool Input::doLength() const
  {return (systemName() == "Isolated") || (systemName() == "Chain") || (systemName() == "Colloid");}
  
  bool Input::doMidFlux() const
  {return dynamicsName() == "BoundaryLangevin";}
  
  bool Input::doSumFlux() const
  {return dynamicsName() == "BoundaryLangevin";}
  
  bool Input::doModiFlux() const
  {return dynamicsName() == "BoundaryLangevin";}
  
  bool Input::doLagrangeMultiplier() const
  {return isConstrained();}
  
  bool Input::doOutThermo() const
  {return true;}  
  
  bool Input::doOutParticles() const
  {return true;}  
  
  bool Input::doXMakeMol() const
  {return systemName() == "NBody" || systemName() == "Bicolor" || systemName() == "Colloid";}    
  
  bool Input::doOutBackUp() const
  {return true;}
  
  bool Input::doOutChain() const
  {return dynamicsName() == "BoundaryLangevin";}
  
  bool Input::doFinalLength() const
  {
    return doLength();
  }  
  
  bool Input::doFinalVelocity() const
  {
    return doVelocity();
  }  
  
  bool Input::doFinalFlux() const
  {
    return dynamicsName() == "BoundaryLangevin"; 
  }  
  
  bool Input::doFinalLagrangeMultiplier() const
  {
    return isConstrained(); 
  }  

}
