#include "LennardJones.hpp"

namespace simol
{
  
  LennardJones::LennardJones(Input const & input):
    Potential(input), 
    epsilon_(input.potentialEpsilon()), 
    sigma_(input.potentialSigma()),
    cutOffRadius_(input.cutOffRatio()*sigma_),
    splineRadius_(input.splineRatio()*cutOffRadius_),
    A_spline_(untruncated(cutOffRadius_)),
    B_spline_(1)
  {}
  
  double LennardJones::splineFunction(double reducedDist) const
  {
     return pow(reducedDist,2)*(A_spline_*reducedDist+B_spline_);
  }

  double LennardJones::untruncated(double dist) const
  { 
    return 4 * epsilon_* (pow(sigma_/dist,12)-pow(sigma_/dist,6));
  }
  
  double LennardJones::operator()(double dist) const
  { 
    return 4 * epsilon_* (pow(sigma_/dist,12)-pow(sigma_/dist,6));
  }
  
  Vector<double> LennardJones::gradient(double dist) const
  { 
    return Vector<double>(1, -24 * epsilon_/ sigma_ * (2*pow(sigma_/dist,13)-pow(sigma_/dist,7)));
  }  
  
  /*Vector<double> LennardJones::gradient(vector<double> const& distVec) const
  {
    double dist = norm(distVec);
    double norm = -24 * epsLJ_/ sigmaLJ_ * (2*pow(sigmaLJ_/dist,13)-pow(sigmaLJ_/dist,7));
    return norm * distVec / dist;
  }*/
  
}