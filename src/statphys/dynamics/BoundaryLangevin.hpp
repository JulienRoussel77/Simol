#ifndef SIMOL_BOUNDARYLANGEVIN_HPP
#define SIMOL_BOUNDARYLANGEVIN_HPP

#include "LangevinBase.hpp"

namespace simol
{
  class BoundaryLangevin : public LangevinBase
  {
  public:
    BoundaryLangevin(Input const& input);
    virtual const double& betaLeft() const;
    virtual const double& betaRight() const;
    virtual const double& temperatureLeft() const;
    virtual const double& temperatureRight() const;
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
    double tauBending_;
  };

}

#endif
