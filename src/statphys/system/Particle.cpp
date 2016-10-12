#ifndef PARTICLE_IMPL_HPP
#define PARTICLE_IMPL_HPP

#include "simol/statphys/system/Particle.hpp"

using std::cout;
using std::endl;

namespace simol
{

  Particle::Particle(double const& mass, DVec const& position0, DVec const& momentum0, double internalEnergy0):
    mass_(mass),
    position_(position0),
    momentum_(momentum0),
    potentialEnergy_(0),
    force_(DVec::Zero(dimension())),
    energyGrad_(DVec::Zero(dimension())),
    energyLapla_(0),
    countdown_(0),
    internalEnergy_(internalEnergy0),
    virial_(0)
  {}

  Particle::Particle(int dimension) :
    Particle(0, DVec::Zero(dimension), DVec::Zero(dimension))
  {}

  Particle::Particle(double const& mass, double const& position0, double const& momentum0, double internalEnergy0):
    Particle(mass, DVec::Constant(1, position0), DVec::Constant(1, momentum0), internalEnergy0)
  {}

  int Particle::dimension() const
  {
    return position_.size();
  }

  double const & Particle::mass() const
  { return mass_; }

  void Particle::resetForce(Potential const& pot)
  {
    potentialEnergy_ = 0;
    force_ = pot.externalForce();
    virial_ = 0;
  }

  //--- primary variables ----

  DVec const & Particle::position() const
  { return position_; }

  DVec & Particle::position()
  { return position_; }

  double const & Particle::position(int i) const
  { return position_(i); }

  double & Particle::position(int i)
  {return position_(i);}

  DVec const & Particle::momentum() const
  { return momentum_; }

  DVec & Particle::momentum()
  { return momentum_; }

  double const & Particle::momentum(int i) const
  { return momentum_(i); }

  double & Particle::momentum(int i)
  { return momentum_(i); }

  //--- functions of the primary variables ----

  double const & Particle::internalEnergy() const
  { return internalEnergy_; }

  double & Particle::internalEnergy()
  { return internalEnergy_; }

  double const & Particle::virial() const
  { return virial_; }

  double & Particle::virial()
  { return virial_; }

  double Particle::kineticEnergy() const
  {
    return pow(momentum_.norm(), 2) / 2 / mass_;
  }

  double const & Particle::potentialEnergy() const
  { return potentialEnergy_; }

  double& Particle::potentialEnergy()
  { return potentialEnergy_; }

  double Particle::energy() const
  { return kineticEnergy() + potentialEnergy_; }

  double Particle::totalEnergyDPDE() const
  { return kineticEnergy() + potentialEnergy_ + internalEnergy_; }

  DVec const& Particle::force() const
  { return force_; }

  DVec& Particle::force()
  { return force_; }

  const double& Particle::force(int i) const
  { return force_(i); }

  double& Particle::force(int i)
  { return force_(i); }

  DVec Particle::velocity() const
  {return momentum_ / mass_;}
  
  double Particle::velocity(int iOfDim) const
  {return momentum_(iOfDim) / mass_;}

  //---- currently specific to chains -----

  DVec const& Particle::energyGrad() const
  { return energyGrad_; }

  DVec& Particle::energyGrad()
  { return energyGrad_; }

  const double& Particle::energyGrad(int i) const
  { return energyGrad_(i); }

  double& Particle::energyGrad(int i)
  { return energyGrad_(i); }

  const double& Particle::energyLapla() const
  {return energyLapla_;}

  double& Particle::energyLapla()
  {return energyLapla_;}

  int const& Particle::countdown() const
  {return countdown_;}

  int& Particle::countdown()
  {return countdown_;}

}

#endif
