#ifndef SIMOL_SPACESINUS_HPP
#define SIMOL_SPACESINUS_HPP

#include "Potential.hpp"

namespace simol
{
  class SpaceSine : public Potential
  {
    public:
      SpaceSine(Input const& input);
      double operator()(DVec const& position) const;
      DVec gradient(DVec const& position) const;
      double laplacian(DVec const& position) const;
      
      double marginalWithoutCoupling(DVec const& position) const;
      double gradientCoupling(DVec const& position) const;

    private:
      double amplitude_;
      double pulsation_;
      double coupling_;
  };

}

#endif