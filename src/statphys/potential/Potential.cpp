#ifndef POTENTIAL_IMPL_HPP
#define POTENTIAL_IMPL_HPP

#include "simol/statphys/potential/Potential.hpp"

#include <cmath>

using std::cout;
using std::endl;

namespace simol
{

  Potential::Potential() {}

  Potential::Potential(Input const & input):
    externalForce_(DVec::Zero(input.dimension()))
  {
    externalForce_(0) = input.externalForce();
    //if (externalForce_(0) != 0)
    // cout << "externalForce = " << externalForce_ << endl;
  }

  ///
  ///Read-only accessor for the external force
  DVec const& Potential::externalForce() const
  {
    return externalForce_;
  }
  ///
  /// Write-read accessor for the external force
  DVec& Potential::externalForce()
  {
    return externalForce_;
  }
  ///
  ///Read-only accessor for the i-th component of the external force
  double const& Potential::externalForce(int const& i) const
  {
    return externalForce_(i);
  }
  ///
  ///Write-read accessor for the i-th component of the external force
  double& Potential::externalForce(int const& i)
  {
    return externalForce_(i);
  }
  
  double const& Potential::parameter1() const
  {
    throw std::runtime_error("parameter1 not defined for this potential");
  }
  
  double const& Potential::parameter2() const
  {
    throw std::runtime_error("parameter2 not defined for this potential");
  }

  double Potential::operator()(DVec const& position) const
  {
    //cout << "Potential::operator()(DVec const& position)" << endl;
    return operator()(position(0));
  }

  double Potential::operator()(double position) const
  {
    //cout << "Potential::operator()(double position)" << endl;
    return operator()(DVec::Constant(1,1,position));
    //return operator()(DVec::Constant(1,1, position));
  }

  double Potential::value(DVec const& position) const
  {
    //cout << "Potential::value(DVec const& position)" << endl;
    return operator()(position);
  }

  double Potential::value(double position) const
  {
    //cout << "Potential::value(double position)" << endl;
    return operator()(position);
  }

  DVec Potential::gradient(DVec const& position) const
  {
    return gradient(position(0));
  }

  DVec Potential::gradient(double position) const
  {
    return -gradient(DVec::Constant(1,1, position));
  }

  DVec Potential::totalForce(DVec const& position) const
  {
    return externalForce_ - gradient(position);
  }

  DVec Potential::totalForce(double position) const
  {
    return externalForce_ - gradient(position);
  }

  DVec Potential::potentialForce(DVec const& position) const
  {
    return - gradient(position);
  }

  DVec Potential::potentialForce(double position) const
  {
    return - gradient(position);
  }

  double Potential::laplacian(DVec const& position) const
  {
    return laplacian(position(0));
  }

  double Potential::laplacian(double position) const
  {
    return laplacian(DVec::Constant(1,1, position));
  }

  ///
  /// ration used in the rejection method
  double Potential::shiftToHarmonic() const
  {
    throw std::invalid_argument("Potential::shiftToHarmonic : Function undefined");
  }
  
  ///
  ///sampling using the rejection method
  //double Potential::drawLaw(double /*localBeta*/, std::shared_ptr<RNG>& /*rng*/) const
  //{
  //  throw std::invalid_argument("Potential::drawLaw : Function undefined");
  //}
  
  ///
  ///sampling using the rejection method
  double Potential::drawLaw(double localBeta, std::shared_ptr<RNG>& rng) const
  {
    double ratio = shiftToHarmonic();
    bool reject = true;
    double xdraw, udraw;
    while (reject)
    {
      xdraw = rng->scalarGaussian() / sqrt(localBeta);
      //cout << ratio << " " << exp(-localBeta * (pow(xdraw, 2)/2 + ratio)) << " " << exp(- localBeta * potential_->value(xdraw)) << endl;
      //cout << xdraw << " " << localBeta * pow(xdraw, 2)/2 - ratio << " >= " << localBeta * potential_->value(xdraw) << endl;
      udraw = rng->scalarUniform();

      reject = (udraw > exp(- localBeta * (value(xdraw) + pow(xdraw, 2) / 2 - ratio)));
      //cout << reject << " " << xdraw << " " << ydraw << endl << endl;
      assert(exp(-localBeta * (pow(xdraw, 2) / 2 - ratio)) >= exp(- localBeta * value(xdraw)));
    }
    return xdraw;
  }


}
#endif
