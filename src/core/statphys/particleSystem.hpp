#ifndef SIMOL_PARTICLESYSTEM_HPP
#define SIMOL_PARTICLESYSTEM_HPP

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
      virtual ~ParticleSystem(){};
      Particle & particle(size_t index);
      std::vector<Particle> & configuration();       
      size_t size() const;
      void launch(Dynamics* model, Output& output, double const& timeStep, int const& numberOfIterations);
      virtual void simulate(Dynamics* model, Output& output, double const& timeStep, int const& numberOfIterations) = 0;
      virtual void computeAllForces(Dynamics const* model) = 0;
      virtual void computeOutput();
      void writeOutput(Output& output);
    protected:
      
      size_t currentTimeIteration_;
      std::vector<Particle> configuration_;

  };

  // redémarrage : nombres aleatoires, juste pannes, tous les N pas de temps
  //
  
  class Isolated : public ParticleSystem
  {
  public:
    Isolated(Input const& input);
    void simulate(Dynamics* model, Output& output, double const& timeStep, int const& numberOfIterations);
    void computeAllForces(Dynamics const* model);
  };
  
    class Chain : public ParticleSystem
  {
  public:
    Chain(Input const& input);
    void simulate(Dynamics* model, Output& output, double const& timeStep, int const& numberOfIterations);
    void computeAllForces(Dynamics const* model);
  };

}






//#include "particleSystem.ipp"

#endif
