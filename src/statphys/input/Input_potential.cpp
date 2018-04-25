#include "simol/statphys/input/Input.hpp"

using std::cout;
using std::endl;
using std::string;

//-- default values for potentials ---
const double defaultPotentialCoeff = 1;
const double defaultInteractionRatio = 1;
const double defaultEpsLJ = 1.;
const double defaultSigma = 1.;
const double defaultSplineRatio = 0.8;
const double defaultCutOffRatio = 2.5;
const double defaultPotentialCenter = 0;
const double defaultCoupling = 0;

namespace simol
{

  string Input::potentialName() const 
  {
    if (pairPotentialName() != "None")
      return pairPotentialName();
    else
      return externalPotentialName();
  }
  
    //return data["Potential"]["Name"].as<string>();}
    
  string Input::externalPotentialName() const
  {
    if (data["Potential"]["External"])
      //if (data["Potential"]["External"]["Name"])
      return data["Potential"]["External"].as<string>();
    else 
      return "None";
  }
  
  string Input::pairPotentialName() const
  {
    if (data["Potential"]["Pair"])
      //if (data["Potential"]["Pair"]["Name"])
      return data["Potential"]["Pair"].as<string>();
    else return "None";
  }
  
  string Input::firstPotentialName() const 
  {
    if (data["Potential"]["FirstType"])
      return data["Potential"]["FirstType"].as<string>();
    else return "None";
  }
  
  string Input::secondPotentialName() const
  {
    if (data["Potential"]["SecondType"])
      return data["Potential"]["SecondType"].as<string>();
    else return "None";
  }
  
  double Input::interactionRatio() const
  {
    if (data["Potential"]["InteractionRatio"])
      return data["Potential"]["InteractionRatio"].as<double>();
    else return defaultInteractionRatio;
  }
  
  bool Input::flag() const
  {
    if (data["Potential"]["Flag"])
      return isYes(data["Potential"]["Flag"].as<string>());
    else return false;
  }
  

  //--- Sinusoidal ---
  double Input::amplitude() const
  {
    if (data["Potential"]["Amplitude"])
      return data["Potential"]["Amplitude"].as<double>();
    else
      return defaultPotentialCoeff;
  }
  
  //--- SpaceSine ---
  
  double Input::coupling() const
  {
    if (data["Potential"]["Coupling"])
      return data["Potential"]["Coupling"].as<double>();
    else
      return defaultCoupling;
  }

  //--- Lennard Jones ---
  double Input::potentialEpsilon() const
  {
    if (data["Potential"]["Epsilon"])
      return data["Potential"]["Epsilon"].as<double>();
    else
      return defaultEpsLJ;
  }

  ///
  /// Represents the equilibrium distance of the potential
  double Input::potentialSigma() const
  {
    if (data["Potential"]["Sigma"])
      return data["Potential"]["Sigma"].as<double>();
    else
      return defaultSigma;
  }
  
  double Input::cutOffRatio() const
  {
    if (data["Potential"]["CutOffRatio"])
      return data["Potential"]["CutOffRatio"].as<double>();
    else
      return defaultCutOffRatio;
  }

  double Input::cutOffRadius() const
  {
    if (data["Potential"]["CutOffRadius"])
      return data["Potential"]["CutOffRadius"].as<double>();
    else
      return cutOffRatio() * potentialSigma();
  }

  double Input::splineRatio() const
  {
    if (data["Potential"]["SplineRatio"])
      return data["Potential"]["SplineRatio"].as<double>();
    else
      return defaultSplineRatio;
  }

  //----- DoubleWell -------
  double Input::height() const
  {
    if (data["Potential"]["Height"])
      return data["Potential"]["Height"].as<double>();
    else throw std::runtime_error("No Height in the input file !");
  }

  double Input::interWell() const
  {
    if (data["Potential"]["InterWell"])
      return data["Potential"]["InterWell"].as<double>();
    else throw std::runtime_error("No InterWell in the input file !");
  }

  //----- Harmonic ----
  double Input::potentialCenter() const
  {
    if (data["Potential"]["Center"])
      return data["Potential"]["Center"].as<double>();
    else
      return defaultPotentialCenter;
  }
  
  double Input::potentialStiffness() const
  {
    if (data["Potential"]["Stiffness"])
      return data["Potential"]["Stiffness"].as<double>();
    else
      return defaultPotentialCoeff;
  }
  
  double Input::potentialDimension() const
  {
    if (data["Potential"]["Dimension"])
      return data["Potential"]["Dimension"].as<double>();
    else
      return dimension();
  }

  //------- FPU ----------
  double Input::potentialAlpha() const
  {
    if (data["Potential"]["Alpha"])
      return data["Potential"]["Alpha"].as<double>();
    else
      return defaultPotentialCoeff;
  }

  ///
  /// The value by default is such that the potential is convex
  double Input::potentialBeta() const
  {
    if (data["Potential"]["Beta"])
      return data["Potential"]["Beta"].as<double>();
    else
      return pow(potentialAlpha(), 2) / (3 * amplitude());
      //return defaultPotentialCoeff;
  }

}
