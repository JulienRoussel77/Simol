#ifndef SIMOL_CHAIN_HPP
#define SIMOL_CHAIN_HPP

#include "simol/statphys/system/System.hpp"
//#include "simol/statphys/dynamics/LangevinBase.hpp"

namespace simol
{

  class Chain : public System
  {
  public:
    Chain(Input const& input);
    //virtual ParticleIterator begin();
    //virtual ParticlePairIterator pairBegin();
    //virtual bool finished(ParticleIterator const& it) const;
    //virtual void incrementeIterator(ParticleIterator& it);    
    
    virtual void incrementePairIterator(ParticlePairIterator& it);
    bool pairFinished(ParticlePairIterator const& it) const;
  };

  class BiChain : public Chain
  {
    public:
      BiChain(Input const& input);
      
      void computeAllForces();
      virtual void computeProfile(Output& output, long int iOfStep) const;
    protected:
      Particle ancorParticle_;
  };

  class TriChain : public Chain
  {
    public:
      TriChain(Input const& input);
      bool const& isOfFixedVolum() const;
      void triInteraction(Particle& particle1, Particle& particle2, Particle& particle3) const;
      void computeAllForces();
      virtual double boundaryPotEnergy() const;
      virtual void computeProfile(Output& output, long int iOfStep) const;
    protected:
      Particle ancorParticle1_;
      Particle ancorParticle2_;
      bool isOfFixedVolum_;
  };

}

#endif
