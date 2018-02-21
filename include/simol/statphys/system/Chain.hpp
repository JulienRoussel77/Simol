#ifndef SIMOL_CHAIN_HPP
#define SIMOL_CHAIN_HPP

#include "simol/statphys/system/System.hpp"
//#include "simol/statphys/dynamics/LangevinBase.hpp"

namespace simol
{

  class Chain : public System
  {
  public:
    Chain(Input const& input, int nbOfWallParticles0 = 0);
    //virtual ParticleIterator begin();
    //virtual ParticlePairIterator pairBegin();
    //virtual bool finished(ParticleIterator const& it) const;
    //virtual void incrementeIterator(ParticleIterator& it);    
    virtual string name() const {return "Chain";}
    
    const Particle& operator()(int iOfParticle = 0) const {return *(configuration_[iOfParticle+nbOfWallParticles_]);}
    Particle& operator()(int iOfParticle = 0) {return *(configuration_[iOfParticle+nbOfWallParticles_]);};

    virtual void incrementePairIterator(ParticlePairIterator& it);
    bool pairFinished(ParticlePairIterator const& it) const;
    virtual void sampleMomenta(DynamicsParameters const& dynaPara);
    
    virtual double length() const;
    virtual double drawPotLaw(double localBeta);
  protected:
    int nbOfWallParticles_;
  };

  class BiChain : public Chain
  {
    public:
      BiChain(Input const& input);
      virtual bool isBiChain() const {return true;}
      
      virtual void samplePositions(DynamicsParameters const& dynaPara);
      void computeAllForces();
      void interaction(Particle& particle1, Particle& particle2) const;
      double leftHeatFlux(int iOfParticle) const;
      double rightHeatFlux(int iOfParticle) const;
      double heatFlux(int iOfParticle) const;
      double heatFluxOnAtom(int iOfParticle) const;            
            
      double computeSumFlux() const;
      void enforceConstraint(double& lagrangeMultiplier, double flux, DynamicsParameters const& dynaPara);
    protected:
  };
  
  class BulkDrivenBiChain : public BiChain
  {
    public:
      BulkDrivenBiChain(Input const& input);
      
      void computeAllForces();
      void interaction(Particle& particle1, Particle& particle2) const;
  };

  class TriChain : public Chain
  {
    public:
      TriChain(Input const& input);
      virtual bool isTriChain() const {return true;}
      
      bool const& isOfFixedVolum() const;
      void triInteraction(Particle& particle1, Particle& particle2, Particle& particle3) const;
      
      virtual void samplePositions(DynamicsParameters const& dynaPara);
      void computeAllForces();
      virtual double boundaryPotEnergy() const;
      int& syst();
    protected:
      //Particle ancorParticle1_;
      //Particle ancorParticle2_;
      bool isOfFixedVolum_;
  };

}

#endif
