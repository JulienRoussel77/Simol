#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <fstream>
#include <vector>

#include "potential.hpp"
//#include "ode/verlet.hpp"
#include "linalg/Vector.hpp"
#include "RNG.hpp"



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
    friend void exact_OU_scheme(Particle & particle, double const gamma, double const beta, double const timeStep, RNG& rng);

    //=============
    // CONSTRUCTORS
    //=============


    public:
      Particle();
      Particle(double const & mass, dvec const & position, dvec const & momentum);
      Particle(double const & mass, double const & position, double const & momentum);

    //==========
    // ACCESSORS
    //==========

    public:

      double const & mass() const;
      dvec const & position() const;
      dvec const & momentum() const;
      double kineticEnergy() const;
      double potentialEnergy(Potential const & potential) const;
      double energy(Potential const & potential) const;

    //=============
    // DATA MEMBERS
    //=============

    private:

      double mass_;
      dvec position_;
      dvec momentum_;
  };

  
  
}



//#include "particle.ipp"





#endif
