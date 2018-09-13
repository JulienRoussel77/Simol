#ifndef SIMOL_LANGEVINBASE_HPP
#define SIMOL_LANGEVINBASE_HPP

#include "Dynamics.hpp"

namespace simol
{
  class LangevinBase : public Dynamics
  {
    public:
      virtual string dynamicsName() const {return "LangevinBase";}
      
      virtual void simulate (System& syst) const; 
      
      virtual const double& gamma() const;
      virtual const double& xi() const;
      int xiNbOfSteps() const;
      virtual bool doMomentaExchange() const;
      virtual void initializeCountdown(Particle& particle) const;
      virtual void updateMomentaExchange(Particle& particle1, Particle& particle2) const;
      virtual void updateOrsteinUhlenbeck(Particle& particle, double localBeta, double localTimeStep) const;
    protected:
      LangevinBase(Input const& input);
      //double gamma_;
      //double xi_;
  };


}


#endif
