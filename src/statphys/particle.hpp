#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "tools.hpp"
#include <fstream>
#include <vector>

#include "potential.hpp"
//#include "ode/verlet.hpp"
#include "linalg/Vector.hpp"
#include "RNG.hpp"



//=====================
// FORWARD DECLARATIONS
//=====================

/*namespace simol
{
  class Particle;
  
    void verlet_scheme(Particle & particle, double timeStep);
    void exact_OU_scheme(Particle & particle, double const gamma, double const beta, double const timeStep, dvec const& randVec);
    void maruyama_scheme(Particle & particle, double const beta_, const double& timeStep, dvec const& randVec);
}*/

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
    /*friend void verlet_scheme(Particle & particle, double timeStep);
    friend void exact_OU_scheme(Particle & particle, double const gamma, double const beta, double const timeStep, dvec const& randVec);
    friend void maruyama_scheme(Particle & particle, double const beta_, const double& timeStep, dvec const& randVec);
*/
    //=============
    // CONSTRUCTORS
    //=============


    public:
      Particle();
      Particle(int dimension);
      Particle(double const & mass, dvec const & position, dvec const & momentum);
      Particle(double const & mass, double const & position, double const & momentum);

    //==========
    // ACCESSORS
    //==========

    public:
      //Particle& operator= (Particle const& particle);
      int dimension() const;
      double const & mass() const;
      dvec const & position() const;
      dvec & position();
      const double& position(int i) const;
      double& position(int i);
      dvec const & momentum() const;
      dvec & momentum();
      const double& momentum(int i) const;
      double& momentum(int i);
      double kineticEnergy() const;
      //double& kineticEnergy();      
      const double& potentialEnergy() const;
      double& potentialEnergy();      
      double energy() const;
      dvec const& force() const;
      dvec& force(); 
      const double& force(size_t i) const;
      double& force(size_t i); 
      dvec const& energyGrad() const;
      dvec& energyGrad();  
      const double& energyGrad(int i) const;
      double& energyGrad(int i);  
			const double& energyLapla() const;
      double& energyLapla();  
      dvec velocity() const;
			int const& countdown() const;
			int& countdown();

    //=============
    // DATA MEMBERS
    //=============

    private:

      double mass_;
      dvec position_;
      dvec momentum_;
      double potentialEnergy_;
      //double kineticEnergy_;
      dvec force_;
      dvec energyGrad_;
			double energyLapla_;
			int countdown_;
  };

  
  
}



//#include "particle.ipp"





#endif
