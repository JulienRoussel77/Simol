#include "simol/statphys/potential/DoubleWell.hpp"

namespace simol
{

  DoubleWell::DoubleWell(Input const & input): 
    Potential(input), 
    height_(input.height()), 
    interWell_(input.interWell()),
    center_(input.potentialCenter())
  {}

  double DoubleWell::operator()(double position) const
  { return height_ * pow(position - center_ - interWell_ / 2, 2) * pow(position - center_ + interWell_ / 2, 2); }

  DVec DoubleWell::gradient(double position) const
  {
    return DVec::Constant(1,1, 4 * height_ * (position - center_) * (position - center_ - interWell_ / 2) * (position - center_ + interWell_ / 2));
  }
  
  double DoubleWell::shiftToHarmonic() const
  {
    assert(height_ > 0);
    return (2 * pow(interWell_, 2)*height_ + 1) / (16 * height_);
  }
  
  DVec DoubleWell::polynomialCoeffs() const
  {
    DVec coeffs(DVec::Zero(5));
    coeffs[0] = pow(center_, 4) - pow(center_*interWell_, 2)/2 + pow(interWell_/2, 4);
    coeffs[1] = -4 *pow(center_, 3)+ center_* pow(interWell_, 2);
    coeffs[2] = 6*pow(center_, 2)- pow(interWell_, 2)/2;
    coeffs[3] = -4 * center_;
    coeffs[4] = 1;
    
    return height_ * coeffs;
  }
}
