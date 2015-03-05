#ifndef PARTICLE_IMPL_HPP
#define PARTICLE_IMPL_HPP

template<class ScalarType> inline
Particle<ScalarType>::Particle(ScalarType const & mass, std::vector<ScalarType> const & positions, std::vector<ScalarType> const & speeds)
:mass_(mass), positions_(positions), speeds_(speeds)
{}

template<class ScalarType> inline
ScalarType const & Particle<ScalarType>::position(size_t const instantIndex) const
{ return positions_[instantIndex]; }

template<class ScalarType> inline
ScalarType const & Particle<ScalarType>::speed(size_t instantIndex) const
{ return speeds_[instantIndex]; }
  
template<class ScalarType> inline
ScalarType const & Particle<ScalarType>::mass() const
{ return mass_; }

#endif
