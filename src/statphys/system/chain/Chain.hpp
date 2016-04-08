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
    //virtual void computeAllForces(Dynamics const& model);
    virtual void thermalize(Dynamics& model);
    virtual void thermalize(LangevinBase& dyna);
    //virtual void computeProfile(Output& output, Dynamics const& model, size_t iOfIteration) const;
    //virtual void writeFinalOutput(Output& output, Dynamics const& model);
  };
  
  class BiChain : public Chain
  {
    Particle ancorParticle_;
  public:
    BiChain(Input const& input);
    void computeAllForces();
    virtual void computeProfile(Output& output, Dynamics const& dyna, size_t iOfIteration) const;
    virtual void computeProfile(Output& output, LangevinBase const& model, size_t iOfIteration) const;
    //virtual void writeFinalOutput(Output& output, Dynamics const& model);
  };
  
  class TriChain : public Chain
  {
    Particle ancorParticle1_;
    Particle ancorParticle2_;
  public:
    TriChain(Input const& input);
    void triInteraction(Particle& particle1, Particle& particle2, Particle& particle3) const;
    void computeAllForces();     
    virtual double boundaryPotEnergy() const;
    virtual void computeProfile(Output& output, Dynamics const& dyna, size_t iOfIteration) const;
    virtual void computeProfile(Output& output, LangevinBase const& model, size_t iOfIteration) const;
    //virtual void writeFinalOutput(Output& output, Dynamics const& model);
  };
  
}

#endif