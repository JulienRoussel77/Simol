#include "simol/statphys/potential/LennardJones.hpp"

namespace simol
{
  //------------- Lennard Jones, truncated and splined to 0 ------
  // define two important distances: cutOffRadius_, after which the potential is set to 0; and splineRadius_, after which the spline acts
  // the spline is sought under the form of a third order polynomial Q(r) = P(x)
  // where x = (r-cutOffRadius_)/(cutOffRadius_-splineRadius_) is the renormalized distance, between -1 and 0, denoted reducedDist in the code
  // and P(x) = x^2 (Ax+B) is such that P(-1) = A_spline_ = LennardJones(splineRadius_)
  //                               and P'(-1) = B_spline_ = LennardJonesDerivative(splineRadius_)*(cutOffRadius_-splineRadius_)
  // a simple computation shows that A*x+B = A_spline_*(3+2*x)+B_spline_*(1+x) [check the values of P(x) and P'(x) at x=-1]

  LennardJones::LennardJones(Input const & input):
    Potential(input),
    epsilon_(input.potentialEpsilon()),
    sigma_(input.potentialSigma()),
    cutOffRadius_(input.cutOffRatio()*sigma_),
    splineRatio_(input.splineRatio()),
    splineRadius_(input.splineRatio()*cutOffRadius_),
    A_spline_(untruncated(splineRadius_)),
    B_spline_(-untruncatedDerivative(splineRadius_) *cutOffRadius_* (1 - splineRatio_)),
    C3_spline_(B_spline_ - 2*A_spline_),
    C4_spline_(3*A_spline_ - B_spline_)
  {}

  //-- functions for spline: Q(r) = P(x) --
  double LennardJones::splineFunction(double dist) const
  {    
    double reducedDist = (1 - dist / cutOffRadius_) / (1 - splineRatio_);
    return pow(reducedDist, 2) * (C3_spline_*reducedDist + C4_spline_);
  }

  double LennardJones::splineFunctionDerivative(double dist) const
  {
    double reducedDist = (1 - dist / cutOffRadius_) / (1 - splineRatio_);
    return -reducedDist / (cutOffRadius_*(1-splineRatio_)) * (3*C3_spline_ * reducedDist + 2*C4_spline_);
  }

  //-- untruncated LennardJones potential --
  double LennardJones::untruncated(double dist) const
  {
    return 4 * epsilon_ * (pow(sigma_ / dist, 12) - pow(sigma_ / dist, 6));
  }

  double LennardJones::untruncatedDerivative(double dist) const
  {
    return -24 * epsilon_ / sigma_ * (2 * pow(sigma_ / dist, 13) - pow(sigma_ / dist, 7));
  }

  //-- actual splined LennardJones: distinguish the domains --
  double LennardJones::operator()(double dist) const
  {
    if (dist < splineRadius_)
      return untruncated(dist);
    else if (dist < cutOffRadius_)
      return splineFunction(dist);
    else return 0;
  }

  DVec LennardJones::gradient(double dist) const
  {
    if (dist < splineRadius_)
      return DVec::Constant(1, untruncatedDerivative(dist));
    else if (dist < cutOffRadius_)
      return DVec::Constant(1, splineFunctionDerivative(dist));
    else return DVec::Constant(1, 0);
  }

  /*DVec LennardJones::gradient(vector<double> const& distVec) const
  {
    double dist = norm(distVec);
    double norm = -24 * epsilon_/ sigma_ * (2*pow(sigma_/dist,13)-pow(sigma_/dist,7));
    return norm * distVec / dist;
  }*/
  
  
  
  LennardJonesRep::LennardJonesRep(Input const & input):
    LennardJones(input)
  {
    cutOffRadius_ = sigma_ * pow(2, 1/6.);
    A_spline_ = untruncated(cutOffRadius_);
  }
  
  double LennardJonesRep::operator()(double dist) const
  {
    if (dist < cutOffRadius_)
      return untruncated(dist) - A_spline_;
    else return 0;
  }

  DVec LennardJonesRep::gradient(double dist) const
  {
    if (dist < cutOffRadius_)
      return DVec::Constant(1, untruncatedDerivative(dist));
    else return DVec::Constant(1, 0);
  }

}
