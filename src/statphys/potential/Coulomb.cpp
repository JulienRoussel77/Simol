#include "simol/statphys/potential/Coulomb.hpp"

namespace simol
{
  //------------- Coulomb potential ------
  
  Coulomb::Coulomb(Input const & input):
    Potential(input),
    epsilon_(input.potentialEpsilon()),
    cutOffRadius_(input.cutOffRatio()),
    coeff_(epsilon_/(1. - 2./sqrt(cutOffRadius_) + 1./cutOffRadius_))
  {
    //cout << "epsilon_ = " << epsilon_ << endl;
    //cout << "coeff_ = " << coeff_ << endl;
  }

  /// -- distinguish the domains --
  ///This potential is not derivable in 0
  double Coulomb::operator()(double dist) const
  {
    if (dist < cutOffRadius_)
      //return coeff_ * (1./dist - 2./sqrt(dist * cutOffRadius_) + 1./cutOffRadius_);
      return coeff_ * pow(1./sqrt(dist) - 1./sqrt(cutOffRadius_), 2);
    else return 0;
  }

  double Coulomb::scalarGradient(double dist) const
  {
    if (dist < cutOffRadius_)
      //return -coeff_ * (pow(dist, 2) - 1./sqrt(pow(dist, 3) * cutOffRadius_));
      return -coeff_/2 * pow(dist, -1.5) * (1./sqrt(dist) - 1./sqrt(cutOffRadius_));
    else return 0;
  }

}
