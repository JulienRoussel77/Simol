#ifndef SIMOL_SYSTEM_HPP
#define SIMOL_SYSTEM_HPP

#include "simol/statphys/Tools.hpp"
#include <vector>
#include "simol/statphys/system/Particle.hpp"
#include "simol/statphys/dynamics/Dynamics.hpp"
#include "simol/statphys/output/Output.hpp"

#include "simol/statphys/potential/AllPotentials.hpp"

#include <iomanip>
using std::setw;

namespace simol
{

  class System
  {
    public:
      System(Input const& input);
      virtual ~System();

      virtual void printName() const;

      //-- fundamental brick: array of particles --
      const std::vector<Particle> & configuration() const;
      std::vector<Particle> & configuration();
      const Particle& getParticle(int index = 0) const;
      Particle& getParticle(int index = 0);
      const int& dimension() const;
      int nbOfParticles() const;

      //-- random numbers --
      const std::shared_ptr<RNG>& rng() const;
      std::shared_ptr<RNG>& rng();

      //-- potential and forces --
      virtual void computeAllForces();
      Potential& potential();
      void interaction(Particle& particle1, Particle& particle2) const;
      double potential(Vector<double> const& position) const;
      double potential(const double& position) const;
      Vector<double> totalForce(Vector<double> const& position) const;
      Vector<double> totalForce(double position) const;
      Vector<double> potentialForce(Vector<double> const & position) const;
      Vector<double> potentialForce(double position) const;
      Vector<double>& externalForce() ;
      Vector<double> const& externalForce() const;
      double& externalForce(const int& i);
      double const& externalForce(const int& i) const;
      double const& potParameter1() const;
      double const& potParameter2() const;

      //-- sampling of initial configurations --
      virtual void thermalize(Dynamics& /*model*/);
      virtual Vector<double> drawMomentum(double localBeta, double mass);
      virtual double drawPotLaw(double localBeta);
      virtual double computeMeanPotLaw(double betaLocal) const;

      //-- output functions --
      virtual void computeProfile(Output& /*output*/, Dynamics const& /*model*/, int /*iOfStep*/)const;

      // currently specific to chains
      virtual double boundaryPotEnergy() const;
      double laplacian(Vector<double> const& position) const;

    protected:
      int dimension_;
      std::vector<Particle> configuration_;
      string settingsPath_;
      std::shared_ptr<RNG> rng_;
      Potential* potential_;
  };



}

#endif
