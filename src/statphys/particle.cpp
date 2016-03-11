#ifndef PARTICLE_IMPL_HPP
#define PARTICLE_IMPL_HPP

#include "particle.hpp"

using std::cout; 
using std::endl; 

namespace simol
{
  //=============
  // CONSTRUCTORS
  //=============
  
   Particle::Particle() :mass_(0), position_(0), momentum_(0), force_(0), energyGrad_(0), countdown_(-1)
   {
    //std::cout << "Particle vide créée !" << std::endl;      
   }
   
   Particle::Particle(int dimension) :
		mass_(0), position_(dimension), 
		momentum_(dimension), force_(dimension), 
		energyGrad_(dimension, 0), energyLapla_(0)
   {
    //std::cout << "Particle vide créée !" << std::endl;      
   }

  Particle::Particle(double const & mass, dvec const & position0, dvec const & momentum0):
		Particle(position0.size())
  {
		mass_ = mass;
		position_ = position0;
		momentum_ = momentum0;
	}
  
  Particle::Particle(double const & mass, double const & position0, double const & momentum0):
    Particle(1)
  {
		mass_ = mass;
		position_ = dvec(1, position0);
		momentum_ = dvec(1, momentum0);
	}
  

  //==========
  // ACCESSORS
  //==========
  
  int Particle::dimension() const
  {
    return position_.size();
  }

  dvec const & Particle::position() const
  { return position_; }
  
  dvec & Particle::position()
  { return position_; }
  
  double const & Particle::position(int i) const
  { return position_(i); }
  
  double & Particle::position(int i)
  { return position_(i); }

  dvec const & Particle::momentum() const
  { return momentum_; }
  
  dvec & Particle::momentum()
  { return momentum_; }
  
  double const & Particle::momentum(int i) const
  { return momentum_(i); }
  
  double & Particle::momentum(int i)
  { return momentum_(i); }

  double const & Particle::mass() const
  { return mass_; }

  double Particle::kineticEnergy() const
  { 
		return pow(momentum_.norm(), 2) / 2 / mass_;
		//return kineticEnergy_; 
	}
  
  /*double& Particle::kineticEnergy()
  { return kineticEnergy_; }*/

  const double& Particle::potentialEnergy() const
  { return potentialEnergy_; }
  
  double& Particle::potentialEnergy()
  { return potentialEnergy_; }
      
  double Particle::energy() const
  { return kineticEnergy() + potentialEnergy_; }
  
  dvec const& Particle::force() const
  { return force_; }
  
  dvec& Particle::force()
  { return force_; }
  
  const double& Particle::force(size_t i) const
  { return force_(i); }
  
  double& Particle::force(size_t i)
  { return force_(i); }
  
  dvec const& Particle::energyGrad() const
  { return energyGrad_; }
  
  dvec& Particle::energyGrad()
  { return energyGrad_; }
  
  const double& Particle::energyGrad(int i) const
  { return energyGrad_(i); }
  
  double& Particle::energyGrad(int i)
  { return energyGrad_(i); }
  
  const double& Particle::energyLapla() const
  {return energyLapla_;}
	
  double& Particle::energyLapla()
	{return energyLapla_;}
  
  dvec Particle::velocity() const
  {return momentum_ / mass_;}
  
  int const& Particle::countdown() const
  {return countdown_;}
	
	int& Particle::countdown()
	{return countdown_;}


}

/*namespace simol
{
  void verlet_scheme(Particle & particle, double timeStep)
  {
    
    particle.momentum_ += timeStep * particle.force_ / 2;
    particle.position_ += timeStep * particle.momentum_ / particle.mass_;
    particle.momentum_ += timeStep * particle.force_ / 2;
  }
  
  void exact_OU_scheme(Particle & particle, double const gamma, double const beta, double const timeStep, dvec const& randVec)
  {
    double alpha = exp(- gamma / particle.mass_ * timeStep);    
    particle.momentum_ = alpha * particle.momentum_ + sqrt((1-pow(alpha, 2))/beta*particle.mass_) * randVec;
  }
  
  void maruyama_scheme(Particle & particle, double const beta_, const double& timeStep, dvec const& randVec)
  {
    particle.position_ += timeStep * particle.force_ + sqrt(2*timeStep/beta_) * randVec;
  }
}*/
#endif
