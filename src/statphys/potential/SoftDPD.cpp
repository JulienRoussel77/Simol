#include "simol/statphys/potential/SoftDPD.hpp"

namespace simol
{
  //------------- soft DPD potential ------
  
  SoftDPD::SoftDPD(Input const & input):
    Potential(input),
    epsilon_(input.potentialEpsilon()),
    cutOffRadius_(input.cutOffRatio())
   {}

  //-- distinguish the domains --
  double SoftDPD::operator()(double dist) const
  {
    if (dist < cutOffRadius_)
      return epsilon_*(1-dist/cutOffRadius_); 
    else return 0;
  }

  DVec SoftDPD::gradient(double dist) const
  {
    if (dist < cutOffRadius_)
      return DVec::Constant(1,1, -epsilon_/cutOffRadius_);
    else return DVec::Constant(1, 0);
  }

}
