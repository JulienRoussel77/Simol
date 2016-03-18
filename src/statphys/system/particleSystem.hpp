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
  
  ParticleSystem* createSystem(Input  const& input);

  class ParticleSystem
  {
    friend ParticleSystem* createSystem(Input  const& input);
    public:
      ParticleSystem(Input const& input);
      virtual ~ParticleSystem();
      const Particle& getParticle(size_t index) const;
      Particle& getParticle(size_t index);
			const size_t& dimension() const;
      const std::vector<Particle> & configuration() const; 
      std::vector<Particle> & configuration();       
      size_t nbOfParticles() const;
      Potential& potential();
      double potential(Vector<double> const& position) const;
      double potential(const double& position) const;
      Vector<double> force(Vector<double> const& position) const;
      double laplacian(Vector<double> const& position) const;
      const std::shared_ptr<RNG> rng() const;
      std::shared_ptr<RNG> rng();
      
      virtual Vector<double> drawMomentum(double localBeta, double mass);
      virtual double drawPotLaw(double localBeta);
      virtual double computeMeanPotLaw(double betaLocal) const;
      
      void launch(Dynamics& model, Output& output);
			virtual void thermalize(Dynamics& /*model*/) {assert(false);};
      virtual void computeAllForces(Dynamics const& model){};
      void computeForce(Particle& particle) const;
      void interaction(Particle& particle1, Particle& particle2) const;
      void triInteraction(Particle& particle1, Particle& particle2, Particle& particle3) const;
			virtual double boundaryPotEnergy() const;
      
			virtual void computeProfile(Output& /*output*/, Dynamics const& /*model*/, size_t /*iOfIteration*/)const{};
      void writeOutput(Output& output, size_t iOfIteration = 0);
      virtual void computeFinalOutput(Output& output, Dynamics const& model);
      virtual void writeFinalOutput(Output& output, Dynamics const& model);
			
    protected:
      size_t dimension_;
      std::vector<Particle> configuration_;
			string settingsPath_;
      std::shared_ptr<RNG> rng_;
      Potential* potential_;
  };

  // redémarrage : nombres aleatoires, juste pannes, tous les N pas de temps
  //
  
  class Isolated : public ParticleSystem
  {
  public:
    Isolated(Input const& input);
    void computeAllForces(Dynamics const& model);
		void writeFinalOutput(Output& output, Dynamics const& model);
  };
	
	class Fluid : public ParticleSystem
  {
  public:
    Fluid(Input const& input);
		void computeAllForces(Dynamics const& model);
		void writeFinalOutput(Output& output, Dynamics const& model);
  };
	
	class Chain : public ParticleSystem
	{
	public:
		Chain(Input const& input);
		virtual void computeAllForces(Dynamics const& model);
		virtual void thermalize(Dynamics& model);
		virtual void computeProfile(Output& output, Dynamics const& model, size_t iOfIteration) const;
		virtual void writeFinalOutput(Output& output, Dynamics const& model);
	};
  
  class BiChain : public Chain
  {
    Particle ancorParticle_;
  public:
    BiChain(Input const& input);
    void computeAllForces(Dynamics const& model);
		virtual void computeProfile(Output& output, Dynamics const& model, size_t iOfIteration) const;
		virtual void writeFinalOutput(Output& output, Dynamics const& model);
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
		virtual void writeFinalOutput(Output& output, Dynamics const& model);
	};

}






//#include "particleSystem.ipp"

#endif
