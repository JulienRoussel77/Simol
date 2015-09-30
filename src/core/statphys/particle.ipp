#ifndef PARTICLE_IMPL_HPP
#define PARTICLE_IMPL_HPP

namespace simol
{
  //=============
  // CONSTRUCTORS
  //=============

  Particle::Particle(double const & mass, 
                                 double const & position,
                                 double const & momentum)
  :mass_(mass), position_(position), momentum_(momentum)
  {}

  //==========
  // ACCESSORS
  //==========

  double const & Particle::position() const
  { return position_; }

  double const & Particle::momentum() const
  { return momentum_; }

  double const & Particle::mass() const
  { return mass_; }

  double Particle::kineticEnergy() const
  { return pow(momentum_, 2)/mass_/2; }

  double Particle::potentialEnergy(Potential const & potential) const
  { return potential(position_); }
      
  double Particle::energy(Potential const & potential) const
  { return kineticEnergy() + potentialEnergy(potential); }

}
#endif
