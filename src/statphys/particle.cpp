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

  Particle::Particle(double const & mass, Vector<double> const & position0, Vector<double> const & momentum0):
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
		position_ = Vector<double>(1, position0);
		momentum_ = Vector<double>(1, momentum0);
	}
  

  //==========
  // ACCESSORS
  //==========
  
  int Particle::dimension() const
  {
    return position_.size();
  }

  Vector<double> const & Particle::position() const
  { return position_; }
  
  Vector<double> & Particle::position()
  { return position_; }
  
  double const & Particle::position(int i) const
  { return position_(i); }
  
  double & Particle::position(int i)
  { return position_(i); }

  Vector<double> const & Particle::momentum() const
  { return momentum_; }
  
  Vector<double> & Particle::momentum()
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
  
  Vector<double> const& Particle::force() const
  { return force_; }
  
  Vector<double>& Particle::force()
  { return force_; }
  
  const double& Particle::force(size_t i) const
  { return force_(i); }
  
  double& Particle::force(size_t i)
  { return force_(i); }
  
  Vector<double> const& Particle::energyGrad() const
  { return energyGrad_; }
  
  Vector<double>& Particle::energyGrad()
  { return energyGrad_; }
  
  const double& Particle::energyGrad(int i) const
  { return energyGrad_(i); }
  
  double& Particle::energyGrad(int i)
  { return energyGrad_(i); }
  
  const double& Particle::energyLapla() const
  {return energyLapla_;}
	
  double& Particle::energyLapla()
	{return energyLapla_;}
  
  Vector<double> Particle::velocity() const
  {return momentum_ / mass_;}
  
  int const& Particle::countdown() const
  {return countdown_;}
	
	int& Particle::countdown()
	{return countdown_;}


}

#endif
