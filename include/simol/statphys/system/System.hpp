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
  
  ///
  ///Implementes an iterator on the pairs of particles
  ///Allows to loop on all the particle pairs easily in every system
  class ParticlePairIterator
  {
  friend class System;
  public:
    ParticlePairIterator(); 
    int const& iOfCell1() const {return iOfCell1_;}
    int& iOfCell1() {return iOfCell1_;}
    Particle const& particle1() const {return **it1_;}
    Particle& particle1() {return **it1_;}
    /*Cell const& cell1() const;
    Cell& cell1();*/
    int const& iOfNeighbor2() const {return iOfNeighbor2_;}
    int& iOfNeighbor2() {return iOfNeighbor2_;}
    /*int const& iOfCell2() const;
    int& iOfCell2();*/
    Particle const& particle2() const {return **it2_;}
    Particle& particle2() {return **it2_;}
    /*Cell const& cell2() const;
    Cell& cell2();*/
    vector<Particle*>::iterator const& endIt2() const {return endIt2_;}
    vector<Particle*>::iterator& endIt2() {return endIt2_;}
    
    //System* syst_;
    vector<Particle*>::iterator it1_;
    //int iOfParticle1_;
    //absolute index of the cell containing particle1
    int iOfCell1_;
    //int iOfParticle2_;
    vector<Particle*>::iterator it2_;
    //contains the neighbor index of the cell containing particle2 relatively to cell1
    int iOfNeighbor2_;
    vector<Particle*>::iterator endIt2_;
  };

  class System
  {
    public:
      System(Input const& input);
      virtual ~System();

      virtual void printName() const;

      //-- fundamental brick: array of particles --
      const std::vector<Particle*> & configuration() const;
      std::vector<Particle*> & configuration();
      const Particle& operator()(int iOfParticle = 0) const {return *(configuration_[iOfParticle]);}
      Particle& operator()(int iOfParticle = 0) {return *(configuration_[iOfParticle]);};
      const Particle& getParticle(int iOfParticle = 0) const;
      Particle& getParticle(int iOfParticle = 0);
      //virtual Particle& getMember(const int& iOfCell, const int& iOfMember) {return getParticle(iOfMember);}
      //virtual Particle const& getMember(const int& iOfCell, const int& iOfMember) const {return getParticle(iOfMember);}
      const int& dimension() const;
      int nbOfParticles() const;
      
      virtual Cell const& cell(int const&) const {throw std::runtime_error("Cell only exist for NBody !");}
      virtual Cell & cell(int const&) {throw std::runtime_error("Cell only exist for NBody !");}

      //-- random numbers --
      const std::shared_ptr<RNG>& rng() const;
      std::shared_ptr<RNG>& rng();
      
      //-- particle iterators --
      /*virtual ParticleIterator begin();
      virtual ConstParticleIterator constBegin();
      virtual ParticleIterator end();*/
      //virtual ConstParticleIterator constBegin();
      /*virtual bool finished(ParticleIteratorBase const& it) const;
      //virtual void incrementeIterator(ParticleIterator& it);  */
      
      //-- particle pair iterators --
      virtual ParticlePairIterator pairBegin();
      virtual bool pairFinished(ParticlePairIterator const& it) const;
      virtual void incrementePairIterator(ParticlePairIterator& it);
      virtual Cell const& getCell2(ParticlePairIterator const&) const {return cell(0);}
      virtual Cell& getCell2(ParticlePairIterator const&) {return cell(0);}

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
      std::vector<Particle*> configuration_;
      //list<Particle> configuration_;
      string settingsPath_;
      std::shared_ptr<RNG> rng_;
      Potential* potential_;
  };



}

#endif
