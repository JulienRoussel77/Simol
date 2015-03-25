#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <fstream>
#include <vector>

#include "Potential.hpp"

#include "core/Vector.hpp"

using namespace simol;


//=====================
// FORWARD DECLARATIONS
//=====================

template<class ScalarType>
class Particle;

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

public:

  //=============
  // CONSTRUCTORS
  //=============

  Particle(ScalarType const & mass, Vector<ScalarType> const & positions, Vector<ScalarType> const & speeds);

  //==========
  // ACCESSORS
  //==========

  ScalarType const & position(size_t instantIndex) const;
  ScalarType const & speed(size_t instantIndex) const;
  ScalarType const & mass() const;

  Vector<ScalarType> const & positions() const
  { return positions_; }

  Vector<ScalarType> const & speeds() const
  { return speeds_; }

private:

  //=============
  // DATA MEMBERS
  //=============

  ScalarType mass_;
  Vector<ScalarType> positions_;
  Vector<ScalarType> speeds_;
};

#include "Particle_impl.hpp"


template<class ScalarType>
void verlet_scheme(Particle<ScalarType> & particle, Potential<ScalarType> const & potential, double delta_t, size_t numberOfIterations)
{
  for (size_t iteration=0; iteration < numberOfIterations; ++iteration)
  {
    ScalarType mass = particle.mass_;
    ScalarType position = particle.positions_[iteration];
    ScalarType speed = particle.speeds_[iteration];

    speed = speed - delta_t * potential.derivative(position) / 2;
    position = position + delta_t * speed / mass;
    speed = speed - delta_t * potential.derivative(position) / 2;

    particle.positions_[iteration+1] = position;
    particle.speeds_[iteration+1] = speed;
  }

}



#endif
