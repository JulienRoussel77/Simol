#ifndef SIMOL_CHAIN_HPP
#define SIMOL_CHAIN_HPP

#include "System.hpp"
#include "LangevinBase.hpp"

namespace simol
{

  class Chain : public System
  {
  public:
    Chain(Input const& input);
    virtual void thermalize(Dynamics& model);
    virtual void thermalize(LangevinBase& dyna);
  };
  
  class BiChain : public Chain
  {
  public:
    BiChain(Input const& input);
    void computeAllForces();
    virtual void computeProfile(Output& output, Dynamics const& dyna, int iOfStep) const;
    virtual void computeProfile(Output& output, LangevinBase const& model, int iOfStep) const;
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
    virtual void computeProfile(Output& output, Dynamics const& dyna, int iOfStep) const;
    virtual void computeProfile(Output& output, LangevinBase const& model, int iOfStep) const;
  protected:
    Particle ancorParticle1_;
    Particle ancorParticle2_;
    bool isOfFixedVolum_;
  };
  
}

#endif