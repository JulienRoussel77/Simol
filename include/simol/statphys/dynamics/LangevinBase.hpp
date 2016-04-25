#ifndef SIMOL_LANGEVINBASE_HPP
#define SIMOL_LANGEVINBASE_HPP

#include "Dynamics.hpp"

namespace simol
{
  class LangevinBase : public Dynamics
  {
    public:
      virtual const double& gamma() const;
      virtual const double& xi() const;
      virtual double& xi();
      int xiNbOfSteps();
      virtual bool doMomentaExchange() const;
      virtual void initializeCountdown(Particle& particle);
      virtual void updateMomentaExchange(Particle& particle1, Particle& particle2);
      void updateOrsteinUhlenbeck(Particle& particle, double localBeta);
    protected:
      LangevinBase(Input const& input);
      double gamma_;
      double xi_;
  };


}


#endif
