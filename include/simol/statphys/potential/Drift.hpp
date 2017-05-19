#ifndef SIMOL_DRIFT_HPP
#define SIMOL_DRIFT_HPP

#include "Potential.hpp"

namespace simol
{

  class Drift : public Potential
  {
    public:
      Drift(Input const& input);
      virtual string classname() const {return "Drift";}
      double operator()(const DVec& position) const;
      DVec gradient(const DVec& position) const;
      double laplacian(double position) const;
      DVec totalForce(const DVec& position, int type) const;
  };

}

#endif