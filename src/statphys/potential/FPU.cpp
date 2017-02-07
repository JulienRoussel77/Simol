#include "simol/statphys/potential/FPU.hpp"

namespace simol
{

  FPU::FPU(Input const & input):
    Potential(input),
    stiffness_(input.potentialStiffness()),
    alpha_(input.potentialAlpha()),
    beta_(input.potentialBeta()),
    qRepartitionFct_(0),
    m1_(0),
    m2_(0),
    m3_(0),
    m4_(0),
    harmonicEquilibrium_(0),
    harmonicStiffness_(0)
  {
    assert(beta_ > 0 || (beta_ == 0 && alpha_ == 0));
    cout << "FPU potential created : V(r) = " << stiffness_ << "r^2/2 + " << alpha_ << "r^3/3 + " << beta_ << "r^4/4" << endl;
    computeHarmonic();
    cout << " m1 = " << m1_ << " m2 = " << m2_ << " m3 = " << m3_ << " m4 = " << m4_ << endl;
    cout << "harmonic potential : Vh(r) = " << harmonicStiffness() << " Dr^2/2 where Delta r = r - r_h and r_h = " << harmonicEquilibrium() << endl;
  }
  
  void FPU::computeHarmonic()
  {
    int nbOfIntegrationNodes = 10000;
    double xMax = 10;
    double step = 2*xMax / nbOfIntegrationNodes;
    for (int iOfNode = 0; iOfNode < nbOfIntegrationNodes; iOfNode++)
    {
      double q = - xMax + iOfNode * step;
      double nuq = exp(-1 * operator()(q)) * step;
      qRepartitionFct_ += nuq;
      m1_ += q * nuq;
      m2_ += pow(q, 2) * nuq;
      m3_ += pow(q, 3) * nuq;
      m4_ += pow(q, 4) * nuq;
    }
    m1_ /= qRepartitionFct_;
    m2_ /= qRepartitionFct_;
    m3_ /= qRepartitionFct_;
    m4_ /= qRepartitionFct_;
    
    double harmonicFirstOrder = (alpha_ * (pow(m2_, 2) - m1_*m3_) + beta_ * (m2_ * m3_ - m1_ * m4_)) / (m2_ - pow(m1_, 2));
    harmonicStiffness_ = stiffness_ + (alpha_ * (m3_ - m1_*m2_) + beta_ * (m4_ - m1_ * m3_)) / (m2_ - pow(m1_, 2));
    harmonicEquilibrium_ = - harmonicFirstOrder / harmonicStiffness_;
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
    //return stiffness_;// + pow(alpha_/stiffness_, 2);
    return harmonicStiffness_;
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
    //return 0;//- alpha_ * stiffness_ / (pow(alpha_, 2) + 3 * pow(stiffness_, 3));
    return harmonicEquilibrium_;
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
