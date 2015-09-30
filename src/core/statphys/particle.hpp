#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <fstream>
#include <vector>

#include "potential.hpp"
//#include "ode/verlet.hpp"
#include "linalg/Vector.hpp"

//=====================
// FORWARD DECLARATIONS
//=====================

namespace simol
{
  class Particle;
}

namespace simol
{

  //==================
  // CLASS DECLARATION
  //==================

  class Particle
  {

    //=================
    // FRIEND FUNCTIONS
    //=================

    //friend void verlet(Particle & particle, HamiltonDynamics const & model, double delta_t);
    friend void verlet_scheme(Particle & particle, Potential const & potential, double timeStep);

    //=============
    // CONSTRUCTORS
    //=============Démarrage : /home/roussel/Travail/simol/build/src/molecular_dynamics -i /home/roussel/Travail/simol/src/test/functional/molecular_dynamics/hamiton.yaml


    public:

      Particle(double const & mass, double const & position, double const & momentum);

    //==========
    // ACCESSORS
    //==========

    public:

      double const & mass() const;
      double const & position() const;
      double const & momentum() const;
      double kineticEnergy() const;
      double potentialEnergy(Potential const & potential) const;
      double energy(Potential const & potential) const;

    //=============
    // DATA MEMBERS
    //=============

    private:

      double mass_;
      double position_;
      double momentum_;
  };

}

#include "particle.ipp"

namespace simol
{
  void verlet_scheme(Particle & particle, Potential const & potential, double timeStep)
  {
    particle.momentum_ -= timeStep * potential.derivative(particle.position_) / 2;
    particle.position_ += timeStep * particle.momentum_ / particle.mass_;
    particle.momentum_ -= timeStep * potential.derivative(particle.position_) / 2;
  }
}



#endif
