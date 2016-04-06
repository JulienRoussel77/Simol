#include "NBody.hpp"

namespace simol
{
  
  // ------------------- NBody ----------------------------
  
NBody::NBody(Input const& input):
    System(input),
    nbOfParticlesPerDimension_(input.nbOfParticlesPerDimension()),
    latticeParameter_(input.latticeParameter()),
    domainSize_(input.latticeParameter()*input.nbOfParticlesPerDimension())
    //DEBUG_(input.outputFolderName()+"debug_potential.txt")
  {
    assert(configuration_.size() > 1);
    for (size_t i = 0; i<input.nbOfParticles(); i++) 
      {
      configuration_[i] = Particle(input.mass(), input.initialPosition(i), input.initialMomentum(i));
      //std::cout << configuration_[i].force() << std::endl;
      }
    // debugage du potentiel !!
    //double pas = 0.005;
    //double Rpot = 2.2; 
    //while (Rpot < 3.1)
    //{
    //  DEBUG_ << Rpot << " " << potential(Rpot) << " " << force(Rpot)(0) << endl;
    //  Rpot += pas;
    //} 
  }
  
  void NBody::printName() const
  {
    cout << "SystemType = NBody" << endl;
  }

  size_t NBody::nbOfParticlesPerDimension() const
  {
    return nbOfParticlesPerDimension_;
  }

  double NBody::latticeParameter() const
  {
    return latticeParameter_;
  }
  
  void NBody::computeAllForces(Dynamics const& dyna)
  {
    //std::cout << "NBody::computeAllForces" << std::endl;
    for (auto&& particle : configuration_)
      dyna.resetForce(particle);
    for (size_t i = 0; i < nbOfParticles(); i++)
      for (size_t j = i+1; j < nbOfParticles(); j++)
      interaction(configuration_[i], configuration_[j]);
  }
    
  void NBody::interaction(Particle& particle1, Particle& particle2) const
  {
    Vector<double> r12 = particle1.position() - particle2.position();
    // take closest periodic image
    double distance = 0.;
    for (int d = 0; d < (int)dimension_; d++)
      {
        r12(d) -= rint(r12(d)/domainSize_)*domainSize_;
        distance += r12(d)*r12(d);
      }
    distance = sqrt(distance);
    // compute energy
    double energy12 = potential(distance);
    // cout << distance << "  " << energy12 << endl;
    particle1.potentialEnergy() += energy12/2;
    particle2.potentialEnergy() += energy12/2;
    // compute forces 
    double force12 = force(distance)(0);
    r12 /= distance;
    particle1.force() += force12 * r12;
    particle2.force() -= force12 * r12;
    // compute pressure, based on the Virial formula for the potential part: P_pot = -\sum_{i < j} r_ij v'(r_ij) / d|Vol|; will divide by d|Vol| at the end
    particle1.virial() += 0.5*force12*distance;
    particle2.virial() += 0.5*force12*distance;
  }
  
}