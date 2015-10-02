#ifndef PARTICLE_IMPL_HPP
#define PARTICLE_IMPL_HPP

#include "particle.hpp"

namespace simol
{
  //=============
  // CONSTRUCTORS
  //=============
  
   Particle::Particle() :mass_(0), position_(0), momentum_(0)
   {
    std::cout << "Particle vide créée !" << std::endl;
      
   }

  Particle::Particle(double const & mass, 
                                 dvec const & position,
                                 dvec const & momentum)
  :mass_(mass), position_(position), momentum_(momentum)
  {}
  
    Particle::Particle(double const & mass, 
                                 double const & position,
                                 double const & momentum)
  :mass_(mass), position_(1), momentum_(1)
  {
    position_(0) = position;
    momentum_(0) = momentum;
  }

  //==========
  // ACCESSORS
  //==========

  dvec const & Particle::position() const
  { return position_; }

  dvec const & Particle::momentum() const
  { return momentum_; }

  double const & Particle::mass() const
  { return mass_; }

  double Particle::kineticEnergy() const
  { return momentum_.norm()/mass_/2; }

  double Particle::potentialEnergy(Potential const & potential) const
  { return potential(position_); }
      
  double Particle::energy(Potential const & potential) const
  { return kineticEnergy() + potentialEnergy(potential); }
  


}

namespace simol
{
  void verlet_scheme(Particle & particle, Potential const & potential, double timeStep)
  {
    /*std::cout << "verlet" << std::endl;
    size_t test = particle.momentum_.size();
        std::cout << "ok" << std::endl;*/
    particle.momentum_ -= timeStep * potential.derivative(particle.position_) / 2;
    particle.position_ += timeStep * particle.momentum_ / particle.mass_;
    particle.momentum_ -= timeStep * potential.derivative(particle.position_) / 2;
  }
  
  void exact_OU_scheme(Particle & particle, double const gamma, double const beta, double const timeStep, dvec const& randVec)
  {
    double alpha = exp(- gamma / particle.mass_ * timeStep);    
    particle.momentum_ = alpha * particle.momentum_ + sqrt((1-pow(alpha, 2))/beta*particle.mass_) * randVec;
  }
}
#endif
