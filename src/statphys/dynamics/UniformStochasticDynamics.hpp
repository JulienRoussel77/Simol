#ifndef SIMOL_UNIFORMSTOCHASTICDYNAMICS_HPP
#define SIMOL_UNIFORMSTOCHASTICDYNAMICS_HPP

#include "StochasticDynamics.hpp"

namespace simol
{
  class UniformStochasticDynamics : public StochasticDynamics
  {
  public:
    UniformStochasticDynamics(Input const& input);
    void initializeMomenta(vector<Particle>& configuration);
    virtual double  temperature() const;
    virtual const double& temperatureLeft() const;
    virtual const double& temperatureRight() const;
    virtual const double& beta() const;
    virtual const double& betaLeft() const;
    virtual const double& betaRight() const;
  protected:
    double beta_;
    double temperature_;
  };

}

#endif
