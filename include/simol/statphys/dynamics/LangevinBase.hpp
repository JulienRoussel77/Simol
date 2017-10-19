#ifndef SIMOL_LANGEVINBASE_HPP
#define SIMOL_LANGEVINBASE_HPP

#include "Dynamics.hpp"

namespace simol
{
  class LangevinBase : public Dynamics
  {
    public:
      virtual string dynamicsName() const {return "LangevinBase";}
      
      virtual const double& gamma() const;
      virtual const double& xi() const;
      int xiNbOfSteps();
      virtual bool doMomentaExchange() const;
      virtual void initializeCountdown(Particle& particle);
      virtual void updateMomentaExchange(Particle& particle1, Particle& particle2);
      virtual void updateOrsteinUhlenbeck(Particle& particle, double localBeta, double localTimeStep);
      virtual void getThermo(Output& output) const;
      virtual void getPressure(Output& output) const;
    protected:
      LangevinBase(Input const& input);
      //double gamma_;
      //double xi_;
  };


}


#endif
