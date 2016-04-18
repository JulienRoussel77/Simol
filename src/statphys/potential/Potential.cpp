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
    externalForce_(input.dimension(), 0)
  {
    externalForce_(0) = input.externalForce();
    //if (externalForce_(0) != 0)
    // cout << "externalForce = " << externalForce_ << endl;
  }

  ///
  ///Read-only accessor for the external force
  Vector<double> const& Potential::externalForce() const
  {
    return externalForce_;
  }
  ///
  /// Write-read accessor for the external force
  Vector<double>& Potential::externalForce()
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

  double Potential::operator()(Vector<double> const& position) const
  {
    //cout << "Potential::operator()(Vector<double> const& position)" << endl;
    return operator()(position(0));
  }

  double Potential::operator()(double position) const
  {
    //cout << "Potential::operator()(double position)" << endl;
    return operator()(Vector<double>(1, position));
  }

  double Potential::value(Vector<double> const& position) const
  {
    //cout << "Potential::value(Vector<double> const& position)" << endl;
    return operator()(position);
  }

  double Potential::value(double position) const
  {
    //cout << "Potential::value(double position)" << endl;
    return operator()(position);
  }

  Vector<double> Potential::gradient(Vector<double> const& position) const
  {
    return gradient(position(0));
  }

  Vector<double> Potential::gradient(double position) const
  {
    return -gradient(Vector<double>(1, position));
  }

  Vector<double> Potential::totalForce(Vector<double> const& position) const
  {
    return externalForce_ - gradient(position);
  }

  Vector<double> Potential::totalForce(double position) const
  {
    return externalForce_ - gradient(position);
  }

  Vector<double> Potential::potentialForce(Vector<double> const& position) const
  {
    return - gradient(position);
  }

  Vector<double> Potential::potentialForce(double position) const
  {
    return - gradient(position);
  }

  double Potential::laplacian(Vector<double> const& position) const
  {
    return laplacian(position(0));
  }

  double Potential::laplacian(double position) const
  {
    return laplacian(Vector<double>(1, position));
  }

  ///
  /// ration used in the rejection method
  double Potential::ratioToHarmonic() const
  {
    throw std::invalid_argument("Potential::ratioToHarmonic : Function undefined");
  }
  ///
  ///sampling using the rejection method
  double Potential::drawLaw(double /*localBeta*/, std::shared_ptr<RNG>& /*rng*/) const
  {
    throw std::invalid_argument("Potential::drawLaw : Function undefined");
  }


}
#endif
