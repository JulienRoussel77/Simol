#include "simol/statphys/potential/Drift.hpp"

namespace simol
{

  Drift::Drift(const Input& input):
    Potential(input)
  {}

  /// This is not a potential force but an exernal drift force, so their is not potential value
  /// The name of the class should probably be changed
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

  ///
  /// the sense of the drift depends on the type of particle (-1, 0 or 1)
  DVec Drift::totalForce(const DVec& /*position*/, int type) const
  {
    return nonEqForce_ * type;
  }

}
