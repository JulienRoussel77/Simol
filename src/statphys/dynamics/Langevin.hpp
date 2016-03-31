#ifndef SIMOL_LANGEVIN_HPP
#define SIMOL_LANGEVIN_HPP

#include "UniformStochasticDynamics.hpp"

namespace simol
{
  class Langevin : public UniformStochasticDynamics
  {
  public:
    Langevin(Input const& input);
    virtual void printName() const;
    virtual const double& gamma() const;
    double sigma() const;
    virtual void updateAfter(Particle& particle);
  protected:
    double gamma_;
  };

}

#endif
