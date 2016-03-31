#ifndef SIMOL_BOUNDARYLANGEVIN_HPP
#define SIMOL_BOUNDARYLANGEVIN_HPP

#include "StochasticDynamics.hpp"

namespace simol
{
  class BoundaryLangevin : public StochasticDynamics
  {
  public:
    BoundaryLangevin(Input const& input);
    virtual const double& betaLeft() const;
    virtual const double& betaRight() const;
		virtual double temperature() const;
    virtual const double& temperatureLeft() const;
    virtual const double& temperatureRight() const;
    double deltaTemperature() const;
		const double& gamma() const;
    double sigmaLeft() const;
    double sigmaRight() const;
    const double& tauBending() const;

    void initializeMomenta(vector<Particle>& configuration);
		virtual void bending(Particle& particle1, Particle& particle2) const;
  protected:
    double betaLeft_;
    double betaRight_;
    double temperatureLeft_;
    double temperatureRight_;
		double gamma_;
    double tauBending_;
  };

}

#endif
