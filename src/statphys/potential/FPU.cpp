#include "simol/statphys/potential/FPU.hpp"

namespace simol
{

  FPU::FPU(Input const & input):
    Potential(input),
    stiffness_(input.potentialStiffness()),
    alpha_(input.potentialAlpha()),
    beta_(input.potentialBeta())
  {
    assert(beta_ > 0 || (beta_ == 0 && alpha_ == 0));
    cout << "FPU potential created : V(r) = " << stiffness_ << "r^2/2 + " << alpha_ << "r^3/3 + " << beta_ << "r^4/4" << endl;
  }


  double FPU::operator()(double distance) const
  { return stiffness_ / 2 * pow(distance, 2) + alpha_ / 3 * pow(distance, 3) + beta_ / 4 * pow(distance, 4); }

  DVec FPU::gradient(double distance) const
  {
    DVec deriv(1);
    deriv(0) = stiffness_ * distance + alpha_ * pow(distance, 2) + beta_ * pow(distance, 3);
    return deriv;
  }

  double FPU::laplacian(double distance) const
  {
    return stiffness_ + 2 * alpha_ * distance + 3 * beta_ * pow(distance, 2);
  }

  double FPU::shiftToHarmonic() const
  {
    assert(stiffness_ == 1);
    if (beta_ > 0)
      return pow(alpha_, 4) / (12 * pow(beta_, 3));
    else
      return 0;
  }
  
  double const& FPU::parameter1() const
  {
    return alpha_;
  }
  
  double const& FPU::parameter2() const
  {
    return beta_;
  }
  
  /// /!\ The mass is supposed to be 1 !
  /// Returns an evaluation of a fitted harmonic force
  double FPU::harmonicForce(double dist) const
  {
    return harmonicStiffness() * (dist - harmonicEquilibrium());
  }
  
  ///
  /// Returns the stiffness of the best fitted harmonic potential
  double FPU::harmonicStiffness() const
  {
    return stiffness_;// + pow(alpha_/stiffness_, 2);
  }
  
  ///
  /// Returns the frequency omega of the best fitted harmonic potential
  double FPU::harmonicFrequency() const
  {
    return sqrt(harmonicStiffness()); //stiffness_ + pow(alpha_, 2);
  }
  
  ///
  /// Returns the equilibrium distance of the best fitted harmonic potential
  double FPU::harmonicEquilibrium() const
  {
    return 0;//- alpha_ * stiffness_ / (pow(alpha_, 2) + 3 * pow(stiffness_, 3));
  }

  /*double FPU::drawLaw(double localBeta, std::shared_ptr<RNG>& rng) const
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

      reject = (udraw > exp(- localBeta * (value(xdraw) + pow(xdraw, 2) / 2 + ratio)));
      //cout << reject << " " << xdraw << " " << ydraw << endl << endl;
      assert(exp(-localBeta * (pow(xdraw, 2) / 2 + ratio)) >= exp(- localBeta * value(xdraw)));
    }
    return xdraw;
  }*/

}
