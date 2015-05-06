#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <fstream>
#include <vector>

#include "Potential.hpp"
#include "ode/verlet.hpp"
#include "linalg/Vector.hpp"

//=====================
// FORWARD DECLARATIONS
//=====================

namespace simol
{
  template<class ScalarType>
  class Particle;
}

namespace simol
{

  //==================
  // CLASS DECLARATION
  //==================

  template<class ScalarType>
  class Particle
  {

    //=================
    // FRIEND FUNCTIONS
    //=================

    friend void verlet<>(Particle<ScalarType> & particle, HamiltonDynamics<ScalarType> const & model, double delta_t);

    //=============
    // CONSTRUCTORS
    //=============

    public:

      Particle(ScalarType const & mass, ScalarType const & position, ScalarType const & momentum);

    //==========
    // ACCESSORS
    //==========

    public:

      ScalarType const & mass() const;
      ScalarType const & position() const;
      ScalarType const & momentum() const;

    //=============
    // DATA MEMBERS
    //=============

    private:

      ScalarType mass_;
      ScalarType position_;
      ScalarType momentum_;
  };

}

#include "Particle.ipp"

namespace simol
{
  template<class ScalarType> inline
  void verlet_scheme(Particle<ScalarType> & particle, Potential<ScalarType> const & potential, double timeStep)
  {
    particle.speed_ -= timeStep * potential.derivative(particle.position_) / 2;
    particle.position_ += timeStep * particle.speed_ / particle.mass_;
    particle.speed_ -= timeStep * potential.derivative(particle.position_) / 2;
  }
}



#endif
