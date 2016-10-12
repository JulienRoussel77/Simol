#include "simol/statphys/input/Input.hpp"

using std::cout;
using std::endl;
using std::string;

const double defaultMass = 1;
const int defaultNbOfParticles = 1;
const int defaultNbOfParticlesPerDimension = 1;
const double defaultLength = 2 * M_PI;
const double defaultLatticeParameter = 1.;
const double defaultHeatCapacity = 1;

namespace simol
{

  string Input::systemName() const
  {
    if (data["System"]["Name"])
      return data["System"]["Name"].as<string>();
    else throw std::runtime_error("No systemName in the input file");
  }

  int Input::dimension() const
  {
    if (data["System"]["Dimension"])
      return data["System"]["Dimension"].as<int>();
    else throw std::runtime_error("No Dimension in the input file !");
  }

  int Input::nbOfParticles() const
  {
    if (data["System"]["Number"])
      return data["System"]["Number"].as<int>();
    else if (data["System"]["NumberPerDimension"])
      return pow(data["System"]["NumberPerDimension"].as<int>(), dimension());
    else return defaultNbOfParticles;
  }

  //--- for NBody systems -----------------
  int Input::nbOfParticlesPerDimension() const
  {
    if (data["System"]["NumberPerDimension"])
      return data["System"]["NumberPerDimension"].as<int>();
    else return defaultNbOfParticlesPerDimension;
  }

  bool Input::doCellMethod() const
  {
    if (data["System"]["Cells"] && data["System"]["Cells"].as<string>() == "yes")
      return true;
    else
      return false;
  }

  double Input::mass() const
  {
    if (data["System"]["Mass"])
      return data["System"]["Mass"].as<double>();
    else return defaultMass;
  }

  double Input::length() const
  {
    if (data["System"]["Length"])
      return data["System"]["Length"].as<double>();
    else
      return defaultLength;
  }

  double Input::latticeParameter() const
  {
    if (data["System"]["LatticeParameter"])
      return data["System"]["LatticeParameter"].as<double>();
    else
      return defaultLatticeParameter;
  }

  // Returns True if the initial conditions are either read from the "settings" file or written in the input file
  bool Input::doSetting() const
  {
    if (data["System"]["Setting"])
      return true;
    return false;
  }
  
  // Returns True if the initial conditions are read from the "settings" file
  bool Input::doRestart() const
  {
    if (data["System"]["Setting"])
      if (data["System"]["Setting"]["RestartFileName"])
        return true;
    return false;
  }



  DVec Input::initialPosition(int const& iOfParticle) const
  {
    DVec q0(DVec::Zero(dimension()));
    if (data["System"]["Position"])
      q0(0) = data["System"]["Position"].as<double>();
    else if (doRestart() && iOfParticle >= 0)
    {
      // the initial positions have been read in the constructor of Input()
      //cout << "using settings for q : " << iOfParticle << "->" << initialPositions_[iOfParticle] << endl;
      q0 = initialPositions_[iOfParticle];
    }
    return q0;
  }

  DVec Input::initialMomentum(int const& iOfParticle) const
  {
    DVec p0(DVec::Zero(dimension()));
    if (data["System"]["Momentum"])
      p0(0) =  data["System"]["Momentum"].as<double>();
    else if (doRestart() && iOfParticle >= 0)
    {
      // the initial positions have been read in the constructor of Input()
      //cout << "using settings for p : " << iOfParticle << "->" << initialMomenta_[iOfParticle] << endl;
      p0 = initialMomenta_[iOfParticle];
    }
    return p0;
  }
  
  double Input::initialInternalEnergy(int const& iOfParticle) const
  {
    if (data["System"]["InternalEnergy"])
      return data["System"]["InternalEnergy"].as<double>();
    else if (doRestart() && iOfParticle >= 0)
    {
      // the initial positions have been read in the constructor of Input()
      //cout << "using settings for e : " << iOfParticle << "->" << initialInternalEnergies_[iOfParticle] << endl;
      return initialInternalEnergies_[iOfParticle];
    }
    else return 1;
  }

  //Chain
  bool Input::isOfFixedVolum() const
  {
    if (data["System"]["Volum"]
        && sameLetters(data["System"]["Volum"].as<string>(), "fixed"))
      return true;

    return false;
  }
  
  //-- heat capacity for DPDE --
  double Input::heatCapacity() const
  {
    if (data["System"]["HeatCapacity"])
      return data["System"]["HeatCapacity"].as<double>();
    else
      return defaultHeatCapacity;
  }

  string Input::restartFileName() const
  {
    if (data["System"]["Setting"])
      if (data["System"]["Setting"]["RestartFileName"])
        return data["System"]["Setting"]["RestartFileName"].as<string>();
    
    return "";
  }

  /*string Input::restartFileName() const
  {
    string name = data["System"]["InitialConditions"].as<string>();
    return name;
  }*/

}
