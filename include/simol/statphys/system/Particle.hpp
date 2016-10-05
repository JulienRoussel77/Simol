#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "simol/statphys/Tools.hpp"
#include "simol/statphys/potential/Potential.hpp"
#include "simol/core/linalg/Vector.hpp"
#include "simol/core/random/RNG.hpp"

namespace simol
{

  class Particle
  {

    public:
      Particle(double const & mass, Vector<double> const & position, Vector<double> const & momentum);
      Particle(int dimension);
      Particle(double const & mass, double const & position, double const & momentum);

      int dimension() const;
      double const & mass() const;
      void resetForce(Potential const& pot);

      //-- access to primary variables --
      Vector<double> const & position() const;
      Vector<double> & position();
      const double& position(int i) const;
      double& position(int i);
      Vector<double> const & momentum() const;
      Vector<double> & momentum();
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
      Vector<double> const& force() const;
      Vector<double>& force();
      const double& force(int i) const;
      double& force(int i);
      double const & virial() const;
      double & virial();
      Vector<double> velocity() const;

      //-- currently specific for chains --
      Vector<double> const& energyGrad() const;
      Vector<double>& energyGrad();
      const double& energyGrad(int i) const;
      double& energyGrad(int i);
      const double& energyLapla() const;
      double& energyLapla();
      int const& countdown() const;
      int& countdown();

    private:

      double mass_;
      Vector<double> position_;
      Vector<double> momentum_;
      double potentialEnergy_;
      Vector<double> force_;
      Vector<double> energyGrad_;
      double energyLapla_;
      int countdown_;
      double internalEnergy_;
      double virial_;   // to compute pressure
  };



}



//#include "particle.ipp"





#endif
