#ifndef SIMOL_PARTICLESYSTEM_HPP
#define SIMOL_PARTICLESYSTEM_HPP

#include "tools.hpp"
#include <vector>
#include "particle.hpp"
#include "dynamics.hpp"
#include "output.hpp"



namespace simol
{
  class ParticleSystem;
  
  ParticleSystem* createSystem(Input  const& input, int const& iOfReplica=0);

  class ParticleSystem
  {
    friend ParticleSystem* createSystem(Input  const& input, int const& iOfReplica);
    public:
      ParticleSystem(Input const& input, int const& /*iOfReplica=0*/);
      virtual ~ParticleSystem(){};
      Particle & getParticle(size_t index);
			const size_t& dimension() const;
      std::vector<Particle> & configuration();       
      size_t nbOfParticles() const;
			virtual void initializeSystem(Dynamics* /*model*/){};
      void launch(Dynamics* model, Output& output);
			virtual void thermalize(Dynamics * /*model*/) {assert(false);};
      virtual void simulate(Dynamics* model);
      virtual void computeAllForces(Dynamics const* model) = 0;
			virtual double boundaryPotEnergy() const;
      virtual void computeOutput(Output& output, Dynamics const* model, size_t iOfIteration);
			virtual void computeProfile(Output& /*output*/, Dynamics const* /*model*/, size_t /*iOfIteration*/){};
      void writeOutput(Output& output, size_t iOfIteration = 0);
      virtual void computeFinalOutput(Output& output, Dynamics const* model);
      virtual void writeFinalOutput(Output& output, Dynamics const* model);
			
    protected:
      size_t dimension_;
      std::vector<Particle> configuration_;
			string settingsPath_;
  };

  // redémarrage : nombres aleatoires, juste pannes, tous les N pas de temps
  //
  
  class Isolated : public ParticleSystem
  {
  public:
    Isolated(Input const& input, int const& iOfReplica=0);
    void computeAllForces(Dynamics const* model);
		void writeFinalOutput(Output& output, Dynamics const* model);
  };
	
	class Fluid : public ParticleSystem
  {
  public:
    Fluid(Input const& input, int const& iOfReplica=0);
		virtual void simulate(Dynamics * model);
		void computeAllForces(Dynamics const* model);
		void writeFinalOutput(Output& output, Dynamics const* model);
  };
	
	class Chain : public ParticleSystem
	{
	public:
		Chain(Input const& input, int const& iOfReplica=0);
		virtual void initializeSystem(Dynamics* model) = 0;
		virtual void computeAllForces(Dynamics const* model) = 0;
		virtual void thermalize(Dynamics * model);
    virtual void simulate(Dynamics * model);
		virtual void computeProfile(Output& output, Dynamics const* model, size_t iOfIteration) = 0;
		virtual void writeFinalOutput(Output& output, Dynamics const* model) = 0;
	};
  
  class BiChain : public Chain
  {
    Particle ancorParticle_;
  public:
    BiChain(Input const& input, int const& iOfReplica=0);
		virtual void initializeSystem(Dynamics* model);
    void computeAllForces(Dynamics const* model);
		virtual void computeProfile(Output& output, Dynamics const* model, size_t iOfIteration);
		virtual void writeFinalOutput(Output& output, Dynamics const* model);
  };
  
  class TriChain : public Chain
  {
    Particle ancorParticle1_;
    Particle ancorParticle2_;
  public:
    TriChain(Input const& input, int const& iOfReplica=0);
		virtual void initializeSystem(Dynamics* model);
    void computeAllForces(Dynamics const* model);		  
		virtual double boundaryPotEnergy() const;
		virtual void computeProfile(Output& output, Dynamics const* model, size_t iOfIteration);
		virtual void writeFinalOutput(Output& output, Dynamics const* model);
	};

}






//#include "particleSystem.ipp"

#endif
