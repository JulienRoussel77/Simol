#ifndef SIMOL_SYSTEM_HPP
#define SIMOL_SYSTEM_HPP

#include "simol/statphys/Tools.hpp"
#include <vector>
#include "simol/statphys/system/Particle.hpp"
#include "simol/statphys/dynamics/Dynamics.hpp"
#include "simol/statphys/output/Output.hpp"

#include "simol/statphys/potential/AllPotentials.hpp"

#include <iomanip>
using std::setw;

namespace simol
{  
  class System;
  class Cell;
  
  class ParticleIterator
  {
  public:
    ParticleIterator(int const& iOfParticle0, Particle* particle0);
    ParticleIterator(System& syst);
    int const& iOfParticle() const {return iOfParticle_;}
    int& iOfParticle() {return iOfParticle_;}
    const Particle& particle() const;
    Particle& particle();
        
    System* syst_;
    int iOfParticle_;
  };
  
  ///
  ///Implementes an iterator on the pairs of particles
  ///Allows to loop on all the particle pairs easily in every system
  class ParticlePairIterator
  {
  public:
    ParticlePairIterator(System& syst);
    int const& iOfParticle1() const {return iOfParticle1_;}
    int& iOfParticle1() {return iOfParticle1_;}
    int const& iOfCell1() const {return iOfCell1_;}
    int& iOfCell1() {return iOfCell1_;}
    Particle const& particle1() const;
    Particle& particle1();
    Cell const& cell1() const;
    Cell& cell1();
    int const& iOfParticle2() const {return iOfParticle2_;}
    int& iOfParticle2() {return iOfParticle2_;}
    int const& iOfNeighbor2() const {return iOfNeighbor2_;}
    int& iOfNeighbor2() {return iOfNeighbor2_;}
    int const& iOfCell2() const;
    int& iOfCell2();
    Particle const& particle2() const;
    Particle& particle2();
    Cell const& cell2() const;
    Cell& cell2();
    
    System* syst_;
    int iOfParticle1_;
    //absolute index of the cell containing particle1
    int iOfCell1_;
    int iOfParticle2_;
    //contains the neighbor index of the cell containing particle2 relatively to cell1
    int iOfNeighbor2_;
  };

  class System
  {
    public:
      System(Input const& input);
      virtual ~System();

      virtual void printName() const;

      //-- fundamental brick: array of particles --
      const std::vector<Particle> & configuration() const;
      std::vector<Particle> & configuration();
      const Particle& getParticle(int index = 0) const;
      Particle& getParticle(int index = 0);
      const int& dimension() const;
      int nbOfParticles() const;
      
      virtual Cell const& cell(int const& /*iOfCell*/) const {throw std::runtime_error("Cell only exist for NBody !");}
      virtual Cell & cell(int const& /*iOfCell*/) {throw std::runtime_error("Cell only exist for NBody !");}

      //-- random numbers --
      const std::shared_ptr<RNG>& rng() const;
      std::shared_ptr<RNG>& rng();
      
      //-- particle iterators --
      virtual ParticleIterator begin();
      virtual bool finished(ParticleIterator const& it) const;
      virtual void incrementeIterator(ParticleIterator& it);  
      
      //-- particle pair iterators --
      virtual ParticlePairIterator pairBegin();
      virtual bool pairFinished(ParticlePairIterator const& it) const;
      virtual void incrementePairIterator(ParticlePairIterator& it);

      //-- potential and forces --
      virtual void computeAllForces();
      Potential& potential();
      void interaction(Particle& particle1, Particle& particle2) const;
      double potential(Vector<double> const& position) const;
      double potential(const double& position) const;
      Vector<double> totalForce(Vector<double> const& position) const;
      Vector<double> totalForce(double position) const;
      Vector<double> potentialForce(Vector<double> const & position) const;
      Vector<double> potentialForce(double position) const;
      Vector<double>& externalForce() ;
      Vector<double> const& externalForce() const;
      double& externalForce(const int& i);
      double const& externalForce(const int& i) const;
      double const& potParameter1() const;
      double const& potParameter2() const;

      //-- sampling of initial configurations --
      virtual Vector<double> drawMomentum(double localBeta, double mass);
      virtual double drawPotLaw(double localBeta);
      virtual double computeMeanPotLaw(double betaLocal) const;

      //-- output functions --
      virtual void computeKineticEnergy(Output& output) const;
      virtual void computePotentialEnergy(Output& output) const;
      virtual void computePressure(Output& output, Dynamics const& dyna) const;
      virtual void computeInternalEnergy(Output& output) const;
      virtual void computeInternalTemperature(Output& output, Dynamics const& dyna) const;
      virtual void computeProfile(Output& /*output*/, Dynamics const& /*model*/, long int /*iOfStep*/)const;

      // currently specific to chains
      virtual double boundaryPotEnergy() const;
      double laplacian(Vector<double> const& position) const;
      

    protected:
      int dimension_;
      std::vector<Particle> configuration_;
      string settingsPath_;
      std::shared_ptr<RNG> rng_;
      Potential* potential_;
  };



}

#endif
