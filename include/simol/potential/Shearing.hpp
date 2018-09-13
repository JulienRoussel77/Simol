#ifndef SIMOL_SPACESINE_HPP
#define SIMOL_SPACESINE_HPP

#include "Potential.hpp"

namespace simol
{
  class Shearing : public Potential
  {
    public:
      Shearing(Input const& input);
      DVec totalForce(DVec const& position) const;
      double operator()(DVec const& position) const;
      DVec gradient(DVec const& position) const;
      double laplacian(DVec const& position) const;

    private:
      double amplitude_;
      double pulsation_;
  };

}

#endif