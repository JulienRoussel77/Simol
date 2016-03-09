#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "tools.hpp"
#include <fstream>
#include <vector>

#include "potential.hpp"
//#include "ode/verlet.hpp"
#include "core/linalg/Vector.hpp"
#include "core/random/RNG.hpp"



//=====================
// FORWARD DECLARATIONS
//=====================

namespace simol
{
  class Particle;

    void verlet_scheme(Particle & particle, double timeStep);
    void exact_OU_scheme(Particle & particle, double const gamma, double const beta, double const timeStep, Vector<double> const& randVec);
    void maruyama_scheme(Particle & particle, double const beta_, double const& timeStep, Vector<double> const& randVec);
}

namespace simol
{

  //==================
  // CLASS DECLARATION
  //==================

  class Particle
  {

    //=================
    // FRIEND FUNCTIONS
    //=================

    //friend void verlet(Particle & particle, HamiltonDynamics const & model, double delta_t);
    friend void verlet_scheme(Particle & particle, double timeStep);
    friend void exact_OU_scheme(Particle & particle, double const gamma, double const beta, double const timeStep, Vector<double> const& randVec);
    friend void maruyama_scheme(Particle & particle, double const beta_, double const& timeStep, Vector<double> const& randVec);

    //=============
    // CONSTRUCTORS
    //=============


    public:
      Particle();
      Particle(int dimension);
      Particle(double const & mass, Vector<double> const & position, Vector<double> const & momentum);
      Particle(double const & mass, double const & position, double const & momentum);

    //==========
    // ACCESSORS
    //==========

    public:
      //Particle& operator= (Particle const& particle);
      int dimension() const;
      double const & mass() const;
      Vector<double> const & position() const;
      Vector<double> & position();
      double const& position(int i) const;
      double& position(int i);
      Vector<double> const & momentum() const;
      Vector<double> & momentum();
      double const& momentum(int i) const;
      double& momentum(int i);
      double const& kineticEnergy() const;
      double& kineticEnergy();
      double const& potentialEnergy() const;
      double& potentialEnergy();
      double energy() const;
      Vector<double> const& force() const;
      Vector<double>& force();
      double const& force(size_t i) const;
      double& force(size_t i);
      Vector<double> const& energyGrad() const;
      Vector<double>& energyGrad();
      double const& energyGrad(int i) const;
      double& energyGrad(int i);
      Vector<double> velocity() const;

    //=============
    // DATA MEMBERS
    //=============

    private:

      double mass_;
      Vector<double> position_;
      Vector<double> momentum_;
      double potentialEnergy_;
      double kineticEnergy_;
      Vector<double> force_;
      Vector<double> energyGrad_;
  };



}



//#include "particle.ipp"





#endif
