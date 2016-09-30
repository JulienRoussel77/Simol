#ifndef SIMOL_LANGEVIN_HPP
#define SIMOL_LANGEVIN_HPP

#include "LangevinBase.hpp"

namespace simol
{
  class Langevin : public LangevinBase
  {
    public:
      Langevin(Input const& input);
      virtual void printName() const;
      double sigma() const;
      //virtual void updateAfter(Particle& particle);
      virtual void computeGeneratorOnBasis(CVBasis& cvBasis, vector<Particle*> const& configuration) const;
    protected:
      double sigma_;
  };

}

#endif
