#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "simol/statphys/Tools.hpp"
#include "simol/statphys/potential/Potential.hpp"
//#include "simol/core/linalg/Vector.hpp"
#include "simol/core/random/RNG.hpp"

namespace simol
{

  class Particle
  {

    public:
      Particle(double const & mass, DVec const & position0, DVec const & momentum0, double internalEnergy0 = 0);
      Particle(int dimension);
      Particle(double const & mass, double const & position0, double const & momentum0, double internalEnergy0 = 0);

      int dimension() const;
      double const & mass() const;
      void resetForce(Potential const& pot);

      //-- access to primary variables --
      DVec const & position() const;
      DVec & position();
      const double& position(int i) const;
      double& position(int i);
      DVec const & momentum() const;
      DVec & momentum();
      const double& momentum(int i) const;
      double& momentum(int i);
      double const & internalEnergy() const;
      double & internalEnergy();
      
      //-- functions depending on the primary variables --
      double kineticEnergy() const;
      const double& potentialEnergy() const;
      double& potentialEnergy();
      double energy() const;
      double totalEnergyDPDE() const;
      DVec const& force() const;
      DVec& force();
      const double& force(int i) const;
      double& force(int i);
      double const & virial() const;
      double & virial();
      DVec velocity() const;
      double velocity(int iOfDim) const;

      //-- currently specific for chains --
      DVec const& energyGrad() const;
      DVec& energyGrad();
      const double& energyGrad(int i) const;
      double& energyGrad(int i);
      const double& energyLapla() const;
      double& energyLapla();
      int const& countdown() const;
      int& countdown();
      
      //-- currently specific for colloids --
      int const & type() const;
      int & type();

    private:

      double mass_;
      DVec position_;
      DVec momentum_;
      double potentialEnergy_;
      DVec force_;
      DVec energyGrad_;
      double energyLapla_;
      int countdown_;
      double internalEnergy_;
      double virial_;   // to compute pressure
      int type_; // for simulations with several populations, like a colloid
  };



}



//#include "particle.ipp"





#endif
