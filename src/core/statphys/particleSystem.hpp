#ifndef SIMOL_PARTICLESYSTEM_HPP
#define SIMOL_PARTICLESYSTEM_HPP

#include <vector>
#include "particle.hpp"
#include "dynamics.hpp"
#include "output.hpp"



namespace simol
{
  class ParticleSystem;
  
  ParticleSystem* createSystem(Input  const& input, int const& indexOfReplica=1);

  class ParticleSystem
  {
    friend ParticleSystem* createSystem(Input  const& input, int const& indexOfReplica);
    public:
      ParticleSystem(Input const& input, int const& indexOfReplica=1);
      virtual ~ParticleSystem(){};
      Particle & particle(size_t index);
      std::vector<Particle> & configuration();       
      size_t numberOfParticles() const;
      void launch(Dynamics* model, Output& output);
      virtual void simulate(Dynamics* model, Output& output);
      virtual void computeAllForces(Dynamics const* model) = 0;
      virtual void computeOutput(Output& output, Dynamics const* model, size_t indexOfIteration);
      void writeOutput(Output& output, size_t indexOfIteration = 0);
      virtual void computeFinalOutput(Output& output, Dynamics const* model);
      void writeFinalOutput(Output& output, Dynamics const* model);
    protected:
      int dimension_;
      std::vector<Particle> configuration_;

  };

  // redémarrage : nombres aleatoires, juste pannes, tous les N pas de temps
  //
  
  class Isolated : public ParticleSystem
  {
  public:
    Isolated(Input const& input, int const& indexOfReplica=1);
    void computeAllForces(Dynamics const* model);
  };
  
    class Chain : public ParticleSystem
  {
  public:
    Chain(Input const& input, int const& indexOfReplica=1);
    void computeAllForces(Dynamics const* model);
    void simulate(Dynamics * model, Output& output);
  };

}






//#include "particleSystem.ipp"

#endif
