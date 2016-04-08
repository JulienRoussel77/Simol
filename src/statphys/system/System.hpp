#ifndef SIMOL_SYSTEM_HPP
#define SIMOL_SYSTEM_HPP

#include "Tools.hpp"
#include <vector>
#include "Particle.hpp"
#include "Dynamics.hpp"
//#include "DPDE.hpp"
#include "Output.hpp"

#include "AllPotentials.hpp"

// fichiers de sortie pour debugage
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
    
    const Particle& getParticle(size_t index = 0) const;
    Particle& getParticle(size_t index = 0);
    const size_t& dimension() const;
    const std::vector<Particle> & configuration() const;
    std::vector<Particle> & configuration();
    size_t nbOfParticles() const;
    Potential& potential();
    double potential(Vector<double> const& position) const;
    double potential(const double& position) const;
    Vector<double> totalForce(Vector<double> const& position) const;
    Vector<double> totalForce(double position) const;
    Vector<double> potentialForce(Vector<double> const & position) const;
    Vector<double> potentialForce(double position) const;
    double laplacian(Vector<double> const& position) const;
    const std::shared_ptr<RNG>& rng() const;
    std::shared_ptr<RNG>& rng();
    
    Vector<double>& externalForce() ;
    Vector<double> const& externalForce() const;
    double& externalForce(const int& i);
    double const& externalForce(const int& i) const;
    
    virtual Vector<double> drawMomentum(double localBeta, double mass);
    virtual double drawPotLaw(double localBeta);
    virtual double computeMeanPotLaw(double betaLocal) const;
    
    //void launch(Dynamics& model, Output& output);
    virtual void thermalize(Dynamics& /*model*/);
    virtual void computeAllForces();
    void interaction(Particle& particle1, Particle& particle2) const;
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



}

#endif
