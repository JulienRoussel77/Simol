#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <fstream>
#include <vector>

#include "Potential.hpp"

#include "LinearAlgebra/Vector.hpp"

//=====================
// FORWARD DECLARATIONS
//=====================

#include "Particle_fwd.hpp"

namespace simol
{

  template<class ScalarType>
  void verlet_scheme(Particle<ScalarType> & particle, Potential<ScalarType> const & potential, double delta_t, size_t numberOfIterations);

  //==================
  // CLASS DECLARATION
  //==================

  template<class ScalarType>
  class Particle
  {

    //=================
    // FRIEND FUNCTIONS
    //=================

    friend void verlet_scheme<>(Particle<ScalarType> & particle, Potential<ScalarType> const & potential, double delta_t, size_t numberOfIterations);

    //=============
    // CONSTRUCTORS
    //=============

    public:

      Particle(ScalarType const & mass, ScalarType const & position, ScalarType const & speed);

    //==========
    // ACCESSORS
    //==========

    public:

      ScalarType const & mass() const;
      ScalarType const & position() const;
      ScalarType const & speed() const;

    //=============
    // DATA MEMBERS
    //=============

    private:

      ScalarType mass_;
      ScalarType position_;
      ScalarType speed_;
  };

}

#include "Particle_impl.hpp"

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
