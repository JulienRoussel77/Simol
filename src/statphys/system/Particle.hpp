#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "Tools.hpp"
#include <fstream>
#include <vector>

#include "Potential.hpp"
//#include "ode/verlet.hpp"
#include "core/linalg/Vector.hpp"
#include "core/random/RNG.hpp"



//=====================
// FORWARD DECLARATIONS
//=====================


namespace simol
{

  //==================
  // CLASS DECLARATION
  //==================

  class Particle
  {

    //=============
    // CONSTRUCTORS
    //=============


    public:
      //Particle();
      Particle(double const & mass, Vector<double> const & position, Vector<double> const & momentum);
      Particle(int dimension);
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
    const double& position(int i) const;
    double& position(int i);
    Vector<double> const & momentum() const;
    Vector<double> & momentum();
    const double& momentum(int i) const;
    double& momentum(int i);
    double const & internalEnergy() const;
    double & internalEnergy();
    double kineticEnergy() const;
    //double& kineticEnergy();
    const double& potentialEnergy() const;
    double& potentialEnergy();
    double energy() const;
    Vector<double> const& force() const;
    Vector<double>& force(); 
    const double& force(size_t i) const;
    double& force(size_t i); 
    Vector<double> const& energyGrad() const;
    Vector<double>& energyGrad();  
    const double& energyGrad(int i) const;
    double& energyGrad(int i);
    const double& energyLapla() const;
    double& energyLapla();  
    Vector<double> velocity() const;
    int const& countdown() const;
    int& countdown();
    
    //=============
    // DATA MEMBERS
    //=============

    private:
    
    double mass_;
    Vector<double> position_;
    Vector<double> momentum_;
    double potentialEnergy_;
    //double kineticEnergy_;
    Vector<double> force_;
    Vector<double> energyGrad_;
    double energyLapla_;
    int countdown_;
    double internalEnergy_;
  };



}



//#include "particle.ipp"





#endif
