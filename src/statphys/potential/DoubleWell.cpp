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

}
