#include "Input.hpp"

using std::cout;
using std::endl;
using std::string;

const double defaultMass = 1;
const int defaultNbOfParticles = 1;
const int defaultNbOfParticlesPerDimension = 1;
const double defaultLength = 2 * M_PI;
const double defaultLatticeParameter = 1.;

namespace simol {

  string Input::systemName() const {return data["System"]["Name"].as<string>();}
  
  int Input::dimension() const 
  { 
    return data["System"]["Dimension"].as<int>();
  }
  
  int Input::nbOfParticles() const {
    if (data["System"]["Number"])
      return data["System"]["Number"].as<int>();
    else if (data["System"]["NumberPerDimension"])
      return pow(data["System"]["NumberPerDimension"].as<int>(),dimension());
    else return defaultNbOfParticles;
  }

  int Input::nbOfParticlesPerDimension() const {
    if (data["System"]["NumberPerDimension"])
      return data["System"]["NumberPerDimension"].as<int>();
    else return defaultNbOfParticlesPerDimension;
  }

  double Input::mass() const {
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
  
  // Returns True if the initial conditions are read from the "settings" file
  bool Input::doFileSettings() const
  {
    if (data["System"]["Settings"])
      if (data["System"]["Settings"].as<string>() == "FileSettings")
	return true;
    return false;
  }
  
  string Input::settingsPath() const
  {
    if (data["System"]["SettingsPath"])
      return simuTypeName()+data["System"]["SettingsPath"].as<string>();
    else return parametersName()+"settings/settings";
  }

  
  Vector<double> Input::initialPosition(int const& iOfParticle) const {
    Vector<double> q0(dimension(), 0);
    if (data["System"]["Position"])
      q0(0) = data["System"]["Position"].as<double>();
    else if (doFileSettings())
      {
	// the initial positions have been read in the constructor of Input()
	cout << "using settings for q : " << iOfParticle << "->" << initialPositions_[iOfParticle] << endl;
	q0 = initialPositions_[iOfParticle];
      }
    return q0;
  }   
  
  Vector<double> Input::initialMomentum(int const& iOfParticle) const 
  {
    Vector<double> p0(dimension(), 0);
    if (data["System"]["Momentum"])
      p0(0) =  data["System"]["Momentum"].as<double>();
    else if (doFileSettings())
      {
	// the initial positions have been read in the constructor of Input()
	cout << "using settings for p : " << iOfParticle << "->" << initialMomenta_[iOfParticle] << endl;
	p0 = initialMomenta_[iOfParticle];
      }
    return p0;
  }  
  
  
}
