#include "simol/statphys/potential/Drift.hpp"

namespace simol
{

  Drift::Drift(const Input& input):
    Potential(input)
  {}


  double Drift::operator()(const DVec& /*position*/) const
  {
    return 0;
  }

  DVec Drift::gradient(const DVec& position) const
  {
    return DVec::Zero(position.rows(), 1);
  }

  double Drift::laplacian(double /*position*/) const
  {
    return 0;
  }

  DVec Drift::totalForce(const DVec& /*position*/, int type) const
  {
    //cout << "Potential::totalForce" << endl;
    //cout << nonEqForce_ << " - " << gradient(position) << " = " << nonEqForce_ - gradient(position) << endl;
    return nonEqForce_ * type;
  }

}
