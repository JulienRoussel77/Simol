#ifndef SIMOL_CHAIN_HPP
#define SIMOL_CHAIN_HPP

#include "System.hpp"


namespace simol
{

  class Chain : public System
  {
  public:
    Chain(Input const& input);
    //virtual void computeAllForces(Dynamics const& model);
    virtual void thermalize(Dynamics& model);
    //virtual void computeProfile(Output& output, Dynamics const& model, size_t iOfIteration) const;
    //virtual void writeFinalOutput(Output& output, Dynamics const& model);
  };
  
  class BiChain : public Chain
  {
    Particle ancorParticle_;
  public:
    BiChain(Input const& input);
    void computeAllForces(Dynamics const& model);
    virtual void computeProfile(Output& output, Dynamics const& model, size_t iOfIteration) const;
    //virtual void writeFinalOutput(Output& output, Dynamics const& model);
  };
  
  class TriChain : public Chain
  {
    Particle ancorParticle1_;
    Particle ancorParticle2_;
  public:
    TriChain(Input const& input);
    void computeAllForces(Dynamics const& model);     
    virtual double boundaryPotEnergy() const;
    virtual void computeProfile(Output& output, Dynamics const& model, size_t iOfIteration) const;
    //virtual void writeFinalOutput(Output& output, Dynamics const& model);
  };
  
}

#endif