#include "FPU.hpp"

namespace simol
{
  
FPU::FPU(Input const & input):
    Potential(input), 
    stiffness_(input.potentialStiffness()),
    alpha_(input.potentialAlpha()),
    beta_(input.potentialBeta())
  {
    assert(beta_ > 0 || (beta_ == 0 && alpha_ == 0));
  }
  
  
  double FPU::operator()(double distance) const
  { return stiffness_/2 * pow(distance, 2) + alpha_/3 * pow(distance, 3) + beta_/4 * pow(distance, 4); }

  Vector<double> FPU::gradient(double distance) const
  { 
    Vector<double> deriv(1);
    deriv(0) = stiffness_ * distance + alpha_ * pow(distance, 2) + beta_ * pow(distance, 3);
    return deriv;
  }
  
  double FPU::laplacian(double distance) const
  { 
    return stiffness_ + 2 * alpha_ * distance + 3 * beta_ * pow(distance, 2);
  }
  
  double FPU::ratioToHarmonic() const
  {
    assert(stiffness_ == 1);
    if (beta_ > 0)
      return - pow(alpha_, 4) / (12 * pow(beta_, 3));
    else
      return 0;
  }
  
  double FPU::drawLaw(double localBeta, std::shared_ptr<RNG>& rng) const
  {
    double ratio = ratioToHarmonic();
    bool reject = true;
    double xdraw, udraw;
    while (reject)
    {
      xdraw = rng->scalarGaussian() / sqrt(localBeta);
      //cout << ratio << " " << exp(-localBeta * (pow(xdraw, 2)/2 + ratio)) << " " << exp(- localBeta * potential_->value(xdraw)) << endl;
      //cout << xdraw << " " << localBeta * pow(xdraw, 2)/2 - ratio << " >= " << localBeta * potential_->value(xdraw) << endl;
      udraw = rng->scalarUniform();
      
      reject = (udraw > exp(- localBeta * (value(xdraw) + pow(xdraw, 2)/2 + ratio)));
      //cout << reject << " " << xdraw << " " << ydraw << endl << endl;
      assert(exp(-localBeta * (pow(xdraw, 2)/2 + ratio)) >= exp(- localBeta * value(xdraw)));
    }
    return xdraw;
  }
  
}