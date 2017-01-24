#ifndef SIMOL_SYSTEM_HPP
#define SIMOL_SYSTEM_HPP

#include "simol/statphys/Tools.hpp"
#include <vector>
#include "simol/statphys/system/Particle.hpp"
//#include "simol/statphys/dynamics/Dynamics.hpp"
//#include "simol/statphys/output/Output.hpp"

#include "simol/statphys/potential/AllPotentials.hpp"
#include "simol/statphys/dynamics/DynamicsParameters.hpp"

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
  
  //class System;
  //System* createSystem(Input const& input);

  class System
  {
  //friend System* createSystem(Input const& input);
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
    //virtual Particle& getMember(const int& iOfCell, const int& iOfMember) {return getParticle(iOfMember);}
    //virtual Particle const& getMember(const int& iOfCell, const int& iOfMember) const {return getParticle(iOfMember);}
    const int& dimension() const;
    int const& nbOfParticles() const;
    
    virtual bool doCells() const {return false;}
    virtual Cell const& cell(int const&) const {throw std::runtime_error("Cell only exist for NBody !");}
    virtual Cell & cell(int const&) {throw std::runtime_error("Cell only exist for NBody !");}
    virtual void reinitializeCells() {}

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
    virtual DVec periodicImage(DVec const& vecDistance) const {return vecDistance;}
    virtual void computeAllForces();
    Potential& potential();
    Potential const& potential() const;
    void interaction(Particle& particle1, Particle& particle2) const;
    double potential(DVec const& position) const;
    double potential(const double& position) const;
    DVec totalForce(DVec const& position) const;
    DVec totalForce(double position) const;
    DVec potentialForce(DVec const & position) const;
    DVec potentialForce(double position) const;
    DVec& externalForce() ;
    DVec const& externalForce() const;
    double& externalForce(const int& i);
    double const& externalForce(const int& i) const;
    double const& potParameter1() const;
    double const& potParameter2() const;

    //-- sampling of initial configurations --
    virtual DVec drawMomentum(double localBeta, double mass);
    virtual double drawPotLaw(double localBeta);
    virtual double computeMeanPotLaw(double betaLocal) const;



    // currently specific to chains
    virtual double boundaryPotEnergy() const;
    double laplacian(DVec const& position) const;
    
    virtual void samplePositions(DynamicsParameters const& dynaPara);
    virtual void sampleMomenta(DynamicsParameters const& dynaPara);

  //protected:
  public:
    int dimension_;
    int nbOfParticles_;
    std::vector<Particle*> configuration_;
    //list<Particle> configuration_;
    string settingsPath_;
    std::shared_ptr<RNG> rng_;
    Potential* potential_;
    bool doSetting_;
  };



}

#endif
