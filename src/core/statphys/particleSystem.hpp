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
    friend ParticleSystem* createSystem(Input  const& input, int const& indexOfReplica);
    public:
      ParticleSystem(Input const& input, int const& /*iOfReplica=0*/);
      virtual ~ParticleSystem(){};
      Particle & getParticle(size_t index);
			const size_t& dimension() const;
      std::vector<Particle> & configuration();       
      size_t numberOfParticles() const;
			virtual void initializeSystem(Dynamics* /*model*/){};
      void launch(Dynamics* model, Output& output);
      virtual void simulate(Dynamics* model);
      virtual void computeAllForces(Dynamics const* model) = 0;
      virtual void computeOutput(Output& output, Dynamics const* model, size_t indexOfIteration);
			virtual void computeProfile(Output& /*output*/, Dynamics const* /*model*/, size_t /*iOfIteration*/){};
      void writeOutput(Output& output, size_t indexOfIteration = 0);
      virtual void computeFinalOutput(Output& output, Dynamics const* model);
      void writeFinalOutput(Output& output, Dynamics const* model);
			
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
  };
  
  class BiChain : public ParticleSystem
  {
    Particle ancorParticle_;
  public:
    BiChain(Input const& input, int const& iOfReplica=0);
    void computeAllForces(Dynamics const* model);
    void simulate(Dynamics * model);
		virtual void computeProfile(Output& output, Dynamics const* model, size_t iOfIteration);
  };
  
  class TriChain : public ParticleSystem
  {
    Particle ancorParticle1_;
    Particle ancorParticle2_;
  public:
    TriChain(Input const& input, int const& iOfReplica=0);
		virtual void initializeSystem(Dynamics* model);
    void computeAllForces(Dynamics const* model);
    void simulate(Dynamics * model);
		virtual void computeProfile(Output& output, Dynamics const* model, size_t iOfIteration);
		void writeFinalOutput(Output& output, Dynamics const* model);
	};

}






//#include "particleSystem.ipp"

#endif
