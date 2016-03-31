#include "StochasticDynamics.hpp"

namespace simol
{
  ///Constructs a purely virtual Dynamics for Dynamics using RNGs
  StochasticDynamics::StochasticDynamics(Input const& input):
    Dynamics(input),
    xi_(input.xi())
  {
		cout << "xi = " << xi() << endl;
	}

  ///Read-only accessor for xi
  const double& StochasticDynamics::xi() const {return xi_;}

  ///Read-write accessor for xi
  double& StochasticDynamics::xi() {return xi_;}
	///
	///Returns the mean number of iterations between 2 random events
	int StochasticDynamics::xiNbOfIterations()
		{return 1 / (xi_ * timeStep_);}
	///
	///Returns true if the dynamics involves a Poisson process (momenta exchange)
	bool StochasticDynamics::doMomentaExchange() const
		{return xi_ > 0;}
	///
	///If the momenta exchange is activated, the times of future events are drawn
	void StochasticDynamics::initializeCountdown(Particle& particle)
	{
		if (doMomentaExchange())
			particle.countdown() = rng_->scalarExponential() * xiNbOfIterations(); // / (xi() * timestep());
		else
			particle.countdown() = -1;
	}
	///
	///Exchanges the momenta of the 2 particles if the time has come
	void StochasticDynamics::updateMomentaExchange(Particle& particle1, Particle& particle2)
	{
		if (particle2.countdown() == 0)
			{
				Vector<double> temp = particle2.momentum();
				particle2.momentum() = particle1.momentum();
				particle1.momentum() = temp;
				particle2.countdown() = rng_->scalarExponential() * xiNbOfIterations(); // / (xi() * timestep());
			}
			else
				particle2.countdown()--;
	}


}
