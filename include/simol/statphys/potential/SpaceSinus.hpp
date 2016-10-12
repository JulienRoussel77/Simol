#ifndef SIMOL_SPACESINUS_HPP
#define SIMOL_SPACESINUS_HPP

#include "Potential.hpp"

namespace simol
{

  class SpaceSinus : public Potential
  {
    public:
      SpaceSinus(Input const& input);
      double operator()(DVec const& position) const;
      DVec gradient(DVec const& position) const;
      double laplacian(DVec const& position) const;

    private:
      double amplitude_;
      double pulsation_;
  };

}

#endif