#ifndef SIMOL_PARTICLESYSTEM_IMPL_HPP
#define SIMOL_PARTICLESYSTEM_IMPL_HPP

#include "System.hpp"

namespace simol
{
  
  System::System(Input const& input):
    dimension_(input.dimension()),
    configuration_(input.nbOfParticles(), Particle(dimension_)),
    settingsPath_(input.settingsPath())
  {
    potential_ = createPotential(input);
  }
  
  ///
  ///Destrucor
  System::~System()
  {
    delete potential_;
  }
  
  void System::printName() const
  {
    cout << "SystemType = System" << endl;
  }

  
  const Particle& System::getParticle(size_t index) const
  { return configuration_[index]; }
  
  Particle& System::getParticle(size_t index) 
  { return configuration_[index]; }
  
  const size_t& System::dimension() const
  {
    return dimension_;
  }
  
  const std::vector< Particle > & System::configuration() const
  { return configuration_; }
  
  std::vector< Particle > & System::configuration() 
  { return configuration_; }
  
  size_t System::nbOfParticles() const
  {
    return configuration_.size();
  }
  
  ///
  ///Returns by value the potential of the dynamics
  Potential& System::potential() {return *potential_;}
  ///
  ///Evaluate the potential for the vector "position"
  double System::potential(Vector<double> const& position) const {return (*potential_)(position);}
  ///
  ///Evaluate the potential for the scalar "position"
  double System::potential(const double& distance) const {return (*potential_)(distance);}
  ///
  ///Evaluate the force for the scalar "position" (potential and external terms)
  Vector<double> System::force(Vector<double> const& position) const
  {
    return potential_->force(position); 
  }
  ///
  ///Evaluate the force for the scalar "position" (potential and external terms)
  Vector<double> System::force(double position) const
  {
    return potential_->force(position); 
  }
  ///
  ///Evaluate the laplacian of the potential for the vector "position"
  double System::laplacian(Vector<double> const& position) const
  {
    return potential_->laplacian(position); 
  }
  
    
  const std::shared_ptr<RNG>& System::rng() const {return rng_;}

  std::shared_ptr<RNG>& System::rng() {return rng_;}
   
  
  double System::boundaryPotEnergy() const
  {return 0;}
  
  ///
  ///Draw a momentum under the invariant measure at inverse temperature "localBeta"
  Vector<double> System::drawMomentum(double localBeta, double mass)
  {
    return sqrt(1 / (localBeta * mass)) * rng_->gaussian();
  }
  ///
  ///Draw a distance or a bending under the invariant measure at inverse temperature "localBeta"
  double System::drawPotLaw(double localBeta)
  {
    return potential_->drawLaw(localBeta, rng_);
  }
  ///Compute the mean distance or bending under the invariant measure
  ///Proceeds to a simple integral quadrature using rectangles
  double System::computeMeanPotLaw(double localBeta) const
  {
    double repFunc = 0;
    double qInteg = 0;
    size_t nbIntegrationNodes = 1000;
    double step = 8. / nbIntegrationNodes;
    Vector<double> deltaQ(1);
    for (size_t iOfNode = 0; iOfNode < nbIntegrationNodes; iOfNode++)
    {
      deltaQ(0) = - 4 + iOfNode * step;
      repFunc += exp(-localBeta * potential(deltaQ));
      qInteg += deltaQ(0) * exp(-localBeta * potential(deltaQ));
    }
    return qInteg / repFunc;
  }
  
  
  ///Computes the force and the energy associated to this pair interaction, and updates these 2 fields
  ///The first 2 derivates of the potential are stored in "particle2"
  void System::interaction(Particle& particle1, Particle& particle2) const
  {
    Vector<double> r12 = particle2.position() - particle1.position();
    double energy12 = potential(r12);
    Vector<double> force12 = force(r12);    // = - v'(q_2 - q_1)
    double lapla12 = laplacian(r12);  // v"(q_2 - q_1)
    
    particle2.potentialEnergy() = energy12;
    particle1.force() -= force12;
    particle2.force() += force12;
    particle2.energyGrad() = -force12;    // v'(q_2 - q_1)
    particle2.energyLapla() = lapla12;    // v"(q_2 - q_1)
  }
  
  ///Computes the force and the energy associated to this triplet interaction, and updates these 2 fields
  ///The first 2 derivates of the potential are stored in "particle2"
  void System::triInteraction(Particle& particle1, Particle& particle2, Particle& particle3) const
  {
    Vector<double> delta = particle3.position() - 2*particle2.position() + particle1.position();
    //double d12 = r12.norm();
    double energy123 = potential(delta);
    Vector<double> force123 = force(delta);    // = - v'(r_2)
    double lapla123 = laplacian(delta);
    
    particle2.potentialEnergy() = energy123;
    particle1.force() += force123;
    particle2.force() -= 2*force123;
    particle3.force() += force123;
    particle2.energyGrad() = -force123;    // - v'(r_2)
    particle2.energyLapla() = lapla123;    // v''(r_2)
  }
  
  
  
  void System::computeProfile(Output& /*output*/, Dynamics const& /*model*/, size_t /*iOfIteration*/) const 
  {throw std::invalid_argument("System::computeProfile : Function undefined");} 
  
  

  
  //void System::computeFinalOutput(Output& /*output*/, Dynamics const& /*dyna*/)
  //{throw std::invalid_argument("System::computeFinalOutput : Function undefined");}
  

  
  
  

  
  
  
  
 
  
}

#endif
