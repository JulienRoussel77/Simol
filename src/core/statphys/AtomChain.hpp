#ifndef SIMOL_ATOMCHAIN_HPP
#define SIMOL_ATOMCHAIN_HPP

//#include "grid.hpp"


namespace simol
{
  template<class ScalarType>
  class AtomChain
  {
    public:
      AtomChain(size_t const numberOfAtoms,
                ScalarType mass,
                ScalarType leftPosition,
                ScalarType leftMomentum,
                ScalarType leftTemperature,
                ScalarType rightTemperature,
                Potential<ScalarType> const & potential);
    public:
      size_t numberOfAtoms() const;
      void simulate(ScalarType const & timeStep);
    private:
      ParticleSystem<ScalarType> particles_;
      Potential<ScalarType> potential_;
      ScalarType leftTemperature_;
      ScalarType rightTemperature_;
  };

  template<class ScalarType> inline
  AtomChain<ScalarType>::AtomChain(size_t const numberOfAtoms, 
                                   ScalarType mass,
                                   ScalarType leftPosition,
                                   ScalarType leftMomentum,
                                   ScalarType leftTemperature,
                                   ScalarType rightTemperature,
                                   Potential<ScalarType> const & potential)
  : particles_(numberOfAtoms, mass, leftPosition, leftMomentum),
    leftTemperature_(leftTemperature),
    rightTemperature_(rightTemperature),
    potential_(potential)
  {}

  template<class ScalarType> inline
  size_t AtomChain<ScalarType>::numberOfAtoms() const
  { return particles_.size(); }
  
  template<class ScalarType>
  void AtomChain<ScalarType>::simulate(ScalarType const & timeStep)
  {
    HamiltonDynamics<ScalarType> model(potential_);
    for(auto & particle : particles_)
    {
      verlet(particle, model, timeStep);
    }
    particles_[1] = particles_[1]; // a completer
    particles_[numberOfAtoms()-1] = particles_[numberOfAtoms()-1]; // a completer
  }


}





#endif
