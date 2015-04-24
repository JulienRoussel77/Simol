#ifndef PARTICLE_IMPL_HPP
#define PARTICLE_IMPL_HPP

namespace simol
{
  //=============
  // CONSTRUCTORS
  //=============

  template<class ScalarType> inline
  Particle<ScalarType>::Particle(ScalarType const & mass, 
                                 ScalarType const & position,
                                 ScalarType const & momentum)
  :mass_(mass), position_(position), momentum_(momentum)
  {}

  //==========
  // ACCESSORS
  //==========

  template<class ScalarType> inline
  ScalarType const & Particle<ScalarType>::position() const
  { return position_; }

  template<class ScalarType> inline
  ScalarType const & Particle<ScalarType>::momentum() const
  { return momentum_; }

  template<class ScalarType> inline
  ScalarType const & Particle<ScalarType>::mass() const
  { return mass_; }

}
#endif
