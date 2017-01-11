#include "simol/statphys/input/Input.hpp"

using std::cout;
using std::endl;
using std::string;

//-- default values for potentials ---
const double defaultPotentialCoeff = 1;
const double defaultEpsLJ = 1.;
const double defaultSigmaLJ = 1.;
const double defaultSplineRatio = 0.8;
const double defaultCutOffRatio = 2.5;

namespace simol
{

  string Input::potentialName() const {return data["Potential"]["Name"].as<string>();}

  //--- Sinusoidal ---
  double Input::amplitude() const
  {
    if (data["Potential"]["Amplitude"])
      return data["Potential"]["Amplitude"].as<double>();
    else
      return defaultPotentialCoeff;
  }

  //--- Lennard Jones ---
  double Input::potentialEpsilon() const
  {
    if (data["Potential"]["Epsilon"])
      return data["Potential"]["Epsilon"].as<double>();
    else
      return defaultEpsLJ;
  }

  double Input::potentialSigma() const
  {
    if (data["Potential"]["Sigma"])
      return data["Potential"]["Sigma"].as<double>();
    else
      return defaultSigmaLJ;
  }

  double Input::cutOffRatio() const
  {
    if (data["Potential"]["CutOffRatio"])
      return data["Potential"]["CutOffRatio"].as<double>();
    else
      return defaultCutOffRatio;
  }

  double Input::splineRatio() const
  {
    if (data["Potential"]["SplineRation"])
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
  double Input::potentialStiffness() const
  {
    if (data["Potential"]["Stiffness"])
      return data["Potential"]["Stiffness"].as<double>();
    else
      return defaultPotentialCoeff;
  }

  //------- FPU ----------
  double Input::potentialAlpha() const
  {
    if (data["Potential"]["Alpha"])
      return data["Potential"]["Alpha"].as<double>();
    else
      return defaultPotentialCoeff;
  }

  double Input::potentialBeta() const
  {
    if (data["Potential"]["Beta"])
      return data["Potential"]["Beta"].as<double>();
    else
      return pow(potentialAlpha(), 2) / (3 * amplitude());
      //return defaultPotentialCoeff;
  }

}
