#include "simol/statphys/system/Colloid.hpp"

namespace simol
{
  //### Colloid ###
  
  Colloid::Colloid(Input const& input, int nbOfColloidParticles0):
    NBody(input),
    nbOfColloidParticles_(nbOfColloidParticles0)
  {
    assert (nbOfColloidParticles_ <= nbOfParticles());
    for (int iOfParticle = nbOfColloidParticles_; iOfParticle < nbOfParticles(); iOfParticle++)
      getParticle(iOfParticle).type() = 1;
    
    //getParticle(45).type() = 1;
    //getParticle(56).type() = 1;
  }
  
  int const& Colloid::nbOfColloidParticles() const
  {return nbOfColloidParticles_;}

  int& Colloid::nbOfColloidParticles()
  {return nbOfColloidParticles_;}
  
  ///Computes the force and the energy associated to this pair interaction, and updates these 2 fields
  ///The first 2 derivates of the potential are stored in "particle2"
  void Colloid::interaction(Particle& particle1, Particle& particle2) const
  {
    int interactionType = ((particle1.type() == 1) && (particle2.type() == 1));
    //if ((particle1.type() == 1) || (particle2.type() == 1))
    //  cout << particle1.type() << " x " << particle2.type() << " -> " << interactionType << endl;
    // take closest periodic image
    DVec r12 = periodicImage(particle1.position() - particle2.position());
    //DVec r12 = particle1.position() - particle2.position();
    double distance = r12.norm();
    // compute energy
    double energy12 = potential(distance, interactionType);

    particle1.potentialEnergy() += energy12 / 2;
    particle2.potentialEnergy() += energy12 / 2;
    // compute forces
    double force12 = potentialForce(distance, interactionType)(0);
    r12 /= distance;
    particle1.force() += force12 * r12;
    particle2.force() -= force12 * r12;
    // compute pressure, based on the Virial formula for the potential part: P_pot = -\sum_{i < j} r_ij v'(r_ij) / d|Vol|; will divide by d|Vol| at the end
    particle1.virial() += 0.5 * force12 * distance;
    particle2.virial() += 0.5 * force12 * distance;
    
    //cout << "inter : " << distance << " " << energy12 << " " << force12 << endl;
    //cout << "--> " << particle1.position().adjoint() << " v " << particle2.position().adjoint() << " -> f= " << particle1.force().adjoint() << endl;
  }
  
  void Colloid::computeAllForces()
  {
    //cout << endl << "------------------computeAllForces----------------------" << endl;
    //for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
    //  cout << "1 : i " << iOfParticle << " " << getParticle(iOfParticle).type() << endl;

    for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
      getParticle(iOfParticle).resetForce(potential());
    
    //for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
    //  cout << "2 : i " << iOfParticle << ", p " << getParticle(iOfParticle).position().adjoint() << ", e " << getParticle(iOfParticle).potentialEnergy() << ", f " << getParticle(iOfParticle).force().adjoint() << endl;
    
    for (ParticlePairIterator it = pairBegin(); !pairFinished(it); incrementePairIterator(it))
      interaction(it.particle1(), it.particle2());
    
    //for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
    //  cout << "3 : i " << iOfParticle << ", p " << getParticle(iOfParticle).position().adjoint() << ", e " << getParticle(iOfParticle).potentialEnergy() << ", f " << getParticle(iOfParticle).force().adjoint() << endl;

  }
  
  ///
  ///Computes the instant value of the observable length
  double Colloid::length() const
  {
    DVec r12 = getParticle(0).position() - getParticle(1).position();
    return r12.norm();
  }
  
}