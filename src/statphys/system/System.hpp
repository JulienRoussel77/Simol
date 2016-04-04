#ifndef SIMOL_PARTICLESYSTEM_HPP
#define SIMOL_PARTICLESYSTEM_HPP

#include "Tools.hpp"
#include <vector>
#include "Particle.hpp"
#include "Dynamics.hpp"
#include "DPDE.hpp"
#include "Output.hpp"



namespace simol
{

  class System
  {
    public:
      System(Input const& input);
      virtual ~System();

      virtual void printName() const;

      const Particle& getParticle(size_t index = 0) const;
      Particle& getParticle(size_t index = 0);
			const size_t& dimension() const;
      const std::vector<Particle> & configuration() const;
      std::vector<Particle> & configuration();
      size_t nbOfParticles() const;
      Potential& potential();
      double potential(Vector<double> const& position) const;
      double potential(const double& position) const;
      Vector<double> force(Vector<double> const& position) const;
      double laplacian(Vector<double> const& position) const;
      const std::shared_ptr<RNG>& rng() const;
      std::shared_ptr<RNG>& rng();

      virtual Vector<double> drawMomentum(double localBeta, double mass);
      virtual double drawPotLaw(double localBeta);
      virtual double computeMeanPotLaw(double betaLocal) const;

      void launch(Dynamics& model, Output& output);
			virtual void thermalize(Dynamics& /*model*/) {assert(false);};
      virtual void computeAllForces(Dynamics const& /*model*/){};
      void interaction(Particle& particle1, Particle& particle2) const;
      void triInteraction(Particle& particle1, Particle& particle2, Particle& particle3) const;
			virtual double boundaryPotEnergy() const;

			virtual void computeProfile(Output& /*output*/, Dynamics const& /*model*/, size_t /*iOfIteration*/)const;
      void writeOutput(Output& output, size_t iOfIteration = 0);
      //virtual void computeFinalOutput(Output& output, Dynamics const& model);
      //virtual void writeFinalOutput(Output& output, Dynamics const& model);

  protected:
    size_t dimension_;
    std::vector<Particle> configuration_;
    string settingsPath_;
    std::shared_ptr<RNG> rng_;
    Potential* potential_;
  };

  class Isolated : public System
  {
  public:
    Isolated(Input const& input);
    void printName() const;
    void computeAllForces(Dynamics const& model);
    //void computeFinalOutput(Output& /*output*/, Dynamics const& /*dyna*/);
    //void writeFinalOutput(Output& output, Dynamics const& model);
  };

  class NBody : public System
  {
  public:
    NBody(Input const& input);
    void computeAllForces(Dynamics const& model);
    //void writeFinalOutput(Output& output, Dynamics const& model);
  };



}






//#include "system.ipp"

#endif
