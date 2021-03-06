#ifndef SIMOL_SYSTEM_HPP
#define SIMOL_SYSTEM_HPP

#include "simol/Tools.hpp"
#include <vector>
#include "simol/system/Particle.hpp"

#include "simol/potential/AllPotentials.hpp"
#include "simol/dynamics/DynamicsParameters.hpp"

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
    int const& iOfNeighbor2() const {return iOfNeighbor2_;}
    int& iOfNeighbor2() {return iOfNeighbor2_;}
    Particle const& particle2() const {return **it2_;}
    Particle& particle2() {return **it2_;}
    vector<Particle*>::iterator const& endIt2() const {return endIt2_;}
    vector<Particle*>::iterator& endIt2() {return endIt2_;}
    
    vector<Particle*>::iterator it1_;
    //absolute index of the cell containing particle1
    int iOfCell1_;
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
    virtual string name() const = 0;
    virtual bool doSetting() const {return doSetting_;}
    virtual bool isBiChain() const {return false;}
    virtual bool isTriChain() const {return false;}
    
    //-- fundamental brick: array of particles --
    const std::vector<Particle*> & configuration() const;
    std::vector<Particle*> & configuration();
    virtual const Particle& operator()(int iOfParticle = 0) const {return *(configuration_[iOfParticle]);}
    virtual Particle& operator()(int iOfParticle = 0) {return *(configuration_[iOfParticle]);};
    const Particle& getParticle(int iOfParticle = 0) const;
    Particle& getParticle(int iOfParticle = 0);
    const int& dimension() const;
    const int& nbOfParticles() const;
    const double& domainSize() const;
    string const& systemSubtype() const;
    
    virtual bool doCells() const {return false;}
    virtual Cell const& cell(int const&) const {throw std::runtime_error("Cell only exist for NBody !");}
    virtual Cell & cell(int const&) {throw std::runtime_error("Cell only exist for NBody !");}
    virtual void reinitializeCells() {}
  private:
    virtual double periodicImage(double position) const;
  public:
    virtual DVec periodicImage(DVec const& vecDistance) const;
    virtual DVec periodicDistance(DVec const& vecDistance) const;

    //-- random numbers --
    const std::shared_ptr<RNG>& rng() const;
    std::shared_ptr<RNG>& rng();
    
    
    //-- particle pair iterators --
    virtual ParticlePairIterator pairBegin();
    virtual bool pairFinished(ParticlePairIterator const& it) const;
    virtual void incrementePairIterator(ParticlePairIterator& it);
    virtual Cell const& getCell2(ParticlePairIterator const&) const {return cell(0);}
    virtual Cell& getCell2(ParticlePairIterator const&) {return cell(0);}

    //-- potential and forces --

    virtual void computeAllForces();
    Potential& externalPotential();
    Potential const& externalPotential() const;
    Potential& pairPotential();
    Potential const& pairPotential() const;
    virtual void interaction(Particle& particle1, Particle& particle2) const;
    double externalPotential(DVec const& position, int type = 0) const;
    double externalPotential(const double& position, int type = 0) const;
    DVec externalForce(DVec const & position, int type = 0) const;
    DVec externalForce(double position, int type = 0) const;
    double const& potParameter1() const;
    double const& potParameter2() const;
    
    DMat forces() const;
    DMat momenta() const;

    //-- sampling of initial configurations --
    virtual DVec drawMomentum(double localBeta, double mass);
    virtual double drawPotLaw(double localBeta);
    virtual double computeMeanPotLaw(double betaLocal) const;



    // currently specific to chains
    virtual double boundaryPotEnergy() const;
    virtual double leftHeatFlux(int /*iOfLink*/) const {throw runtime_error("leftHeatFlux not implemented for this system !");};
    virtual double rightHeatFlux(int /*iOfLink*/) const {throw runtime_error("rightHeatFlux not implemented for this system !");};
    virtual double heatFlux(int /*iOfLink*/) const {throw runtime_error("heatFlux not implemented for this system !");};
    virtual double heatFluxOnAtom(int /*iOfParticle*/) const {throw runtime_error("heatFluxOnAtom not implemented for this system !");};
    
    // currently specific to bicolor systems
    virtual void enforceConstraint(double /*fixedVelocity*/, DynamicsParameters const& /*dynaPara*/, bool /*updateLagrangeMultiplier*/) {throw runtime_error("enforceConstraint not implemented for this system !");}
    
    virtual void samplePositions(DynamicsParameters const& dynaPara);
    virtual void sampleMomenta(DynamicsParameters const& dynaPara);
    virtual void sampleInternalEnergies();
    
    virtual double length() const;
    virtual double velocity() const;
    virtual double force() const;
    
    virtual double& lagrangeMultiplier() {return lagrangeMultiplier_;}
    virtual const double& lagrangeMultiplier() const {return lagrangeMultiplier_;}
  protected:
    int dimension_;
    // nomber of particles
    int nbOfParticles_;
    // vector containing pointers to each particle
    std::vector<Particle*> configuration_;
    //  ugly: contains a path to the address of the file giving the inital configuration (if it exists)
    string settingsPath_;
    // pointer to the random number generator
    std::shared_ptr<RNG> rng_;
    // part of the potential which is not an interaction potential
    // eg a confining potential or an external drift force
    Potential* externalPotential_;
    // pair interaction potential (eg LJ)
    Potential* pairPotential_;
    // Is true if the system is initialized using a file
    bool doSetting_;
    // size of the domain, can be infinity
    double domainSize_;
    // used for the Bicolor system and indicates the type of drift (OneDrift, TwoDrift or ColorDrift)
    string systemSubtype_;
    // for constrained systems, instantaneous value of the Lagrange multiplier
    double lagrangeMultiplier_;
  };



}

#endif
