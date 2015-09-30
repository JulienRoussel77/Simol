#ifndef SIMOL_ATOMCHAIN_HPP
#define SIMOL_ATOMCHAIN_HPP

//#include "grid.hpp"


namespace simol
{
  class AtomChain
  {
    public:
      AtomChain(size_t const numberOfAtoms,
                double mass,
                double leftPosition,
                double leftMomentum,
                double leftTemperature,
                double rightTemperature,
                Potential const & potential);
    public:
      size_t numberOfAtoms() const;
      void simulate(double const & timeStep);
    private:
      ParticleSystem particles_;
      Potential potential_;
      double leftTemperature_;
      double rightTemperature_;
  };

  AtomChain::AtomChain(size_t const numberOfAtoms, 
                                   double mass,
                                   double leftPosition,
                                   double leftMomentum,
                                   double leftTemperature,
                                   double rightTemperature,
                                   Potential const & potential)
  : particles_(numberOfAtoms, mass, leftPosition, leftMomentum),
    leftTemperature_(leftTemperature),
    rightTemperature_(rightTemperature),
    potential_(potential)
  {}

  size_t AtomChain::numberOfAtoms() const
  { return particles_.size(); }
  
  /*void AtomChain::simulate(double const & timeStep)
  {
    /*Dynamics* model = createDynamics(potential_);
    for(auto & particle : particles_.configuration())
    {
      verlet_scheme(particle, model->potential(), timeStep);
    }
    //particles_[1] = particles_[1]; // a completer
    //particles_[numberOfAtoms()-1] = particles_[numberOfAtoms()-1]; // a completer
  }*/


}





#endif
