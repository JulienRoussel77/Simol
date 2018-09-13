#include "simol/potential/DoubleWell.hpp"

namespace simol
{

  DoubleWell::DoubleWell(Input const & input): 
    Potential(input), 
    interWell_(input.interWell()),
    coeff_(input.height() * pow(2/interWell_, 4))
  {}

  double DoubleWell::operator()(double position) const
  { 
    //cout << coeff_ << " " << interWell_ << " " << center_ << endl;
    //ofstream test("testest", std::ofstream::app);
    //test << position << " " << coeff_ * pow(position - center_ - interWell_ / 2, 2) * pow(position - center_ + interWell_ / 2, 2) << endl;
    return coeff_ * pow(position - center_ - interWell_ / 2, 2) * pow(position - center_ + interWell_ / 2, 2);
  }
  
  double DoubleWell::symmetricValue(double position) const
  {
    return value(position);
  }
  
  double DoubleWell::skewsymmetricValue(double /*position*/) const
  {
    return 0;
  }

  double DoubleWell::scalarGradient(double position) const
  {
    return 4 * coeff_ * (position - center_) * (position - center_ - interWell_ / 2) * (position - center_ + interWell_ / 2);
  }
  
  double DoubleWell::laplacian(double position) const
  {
    return 4 * coeff_ * (3*pow(position-center_, 2) - pow(interWell_ / 2, 2));
  }
  
  double DoubleWell::shiftToHarmonic() const
  {
    assert(coeff_ > 0);
    return (2 * pow(interWell_, 2)*coeff_ + 1) / (16 * coeff_);
  }
  
  DVec DoubleWell::polyCoeffs() const
  {
    DVec coeffs(DVec::Zero(5));
    /*coeffs[0] = pow(center_, 4) - pow(center_*interWell_, 2)/2 + pow(interWell_/2, 4);
    coeffs[1] = -4 *pow(center_, 3)+ center_* pow(interWell_, 2);
    coeffs[2] = 6*pow(center_, 2)- pow(interWell_, 2)/2;
    coeffs[3] = -4 * center_;
    coeffs[4] = 1;*/
    
    coeffs[0] = pow(interWell_/2, 4);
    coeffs[2] = -pow(interWell_, 2)/2;
    coeffs[4] = 1;
    
    return coeff_ * coeffs;
  }
}
