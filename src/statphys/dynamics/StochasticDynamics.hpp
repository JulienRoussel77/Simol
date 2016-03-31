#ifndef STOCHASTICDYNAMICS_HPP
#define STOCHASTICDYNAMICS_HPP

#include "Dynamics.hpp"

namespace simol
{
  class StochasticDynamics : public Dynamics
  {
    public:
        StochasticDynamics(Input const& input);
	    virtual const double& xi() const;
	    virtual double& xi();
	    int xiNbOfIterations();
	    virtual bool doMomentaExchange() const;
	    virtual void initializeCountdown(Particle& particle);
	    virtual void updateMomentaExchange(Particle& particle1, Particle& particle2);
    private:
	    double xi_;
  };


}


#endif
