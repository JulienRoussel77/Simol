#ifndef SIMOL_CHAIN_HPP
#define SIMOL_CHAIN_HPP

#include "simol/system/System.hpp"
//#include "simol/dynamics/LangevinBase.hpp"

namespace simol
{


  class Chain : public System
  {
    public:
      Chain(Input const& input);
      
      virtual string name() const {return "Chain";}
    
      const Particle& operator()(int iOfParticle = 0) const {return *(configuration_[iOfParticle]);}
      Particle& operator()(int iOfParticle = 0) {return *(configuration_[iOfParticle]);};

      virtual void incrementePairIterator(ParticlePairIterator& it);
      bool pairFinished(ParticlePairIterator const& it) const;
      virtual void sampleMomenta(DynamicsParameters const& dynaPara);
      
      virtual double length() const;
      virtual double drawPotLaw(double localBeta);
      virtual bool isBiChain() const {return true;}
      
      virtual void samplePositions(DynamicsParameters const& dynaPara);
      void computeAllForces();
      void interaction(Particle& particle1, Particle& particle2) const;
      double leftHeatFlux(int iOfParticle) const;
      double rightHeatFlux(int iOfParticle) const;
      double heatFlux(int iOfParticle) const;
      double heatFluxOnAtom(int iOfParticle) const;            
            
      double computeSumFlux() const;
      void enforceConstraint(double flux, DynamicsParameters const& dynaPara, bool updateLagrangeMultiplier);
    protected:
  };
  
  class BulkDrivenChain : public Chain
  {
    public:
      BulkDrivenChain(Input const& input);
      
      void computeAllForces();
      void interaction(Particle& particle1, Particle& particle2) const;
  };


}

#endif
