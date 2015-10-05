#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <fstream>
#include <vector>

#include "potential.hpp"
//#include "ode/verlet.hpp"
#include "linalg/Vector.hpp"
#include "RNG.hpp"



//=====================
// FORWARD DECLARATIONS
//=====================

namespace simol
{
  class Particle;
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
    friend void verlet_scheme(Particle & particle, Potential const & potential, double timeStep);
    friend void exact_OU_scheme(Particle & particle, double const gamma, double const beta, double const timeStep, dvec const& randVec);

    //=============
    // CONSTRUCTORS
    //=============


    public:
      Particle();
      Particle(double const & mass, dvec const & position, dvec const & momentum);
      Particle(double const & mass, double const & position, double const & momentum);

    //==========
    // ACCESSORS
    //==========

    public:

      double const & mass() const;
      dvec const & position() const;
      dvec const & momentum() const;
      double const& kineticEnergy() const;
      double& kineticEnergy();      
      double const& potentialEnergy() const;
      double& potentialEnergy();      
      double const& energy() const;
      dvec const& force() const;
      dvec& force();  

    //=============
    // DATA MEMBERS
    //=============

    private:

      double mass_;
      dvec position_;
      dvec momentum_;
      double potentialEnergy_;
      double kineticEnergy_;
      dvec force_;
  };

  
  
}



//#include "particle.ipp"





#endif
