#include "simol/statphys/system/Colloid.hpp"

namespace simol
{
  //### Colloid ###
  
  Colloid::Colloid(Input const& input, int nbOfColloidParticles0):
    NBody(input),
    nbOfColloidParticles_(nbOfColloidParticles0)
  {
    assert (nbOfColloidParticles_ <= nbOfParticles());
    for (int iOfParticle = 0; iOfParticle < nbOfColloidParticles(); iOfParticle++)
      getParticle(iOfParticle).type() = 1;
    
  }
  
  int const& Colloid::nbOfColloidParticles() const
  {return nbOfColloidParticles_;}

  int& Colloid::nbOfColloidParticles()
  {return nbOfColloidParticles_;}
  
  ///Computes the force and the energy associated to this pair interaction, and updates these 2 fields
  ///The first 2 derivates of the potential are stored in "particle2"
  void Colloid::interaction(Particle& particle1, Particle& particle2) const
  {
    
    int interactionType = particle1.type() + particle2.type();
    if (interactionType == 2) return;
    
    // take closest periodic image
    DVec r12 = periodicDistance(particle2.position() - particle1.position());
    
    double distance = r12.norm();
    r12 /= distance;
    // compute energy
    double energy12 = pairPotential()(distance, interactionType);

    particle1.potentialEnergy() += energy12 / 2;
    particle2.potentialEnergy() += energy12 / 2;
    // compute forces (positive -> repulsive / negative -> attractive)
    double force12 = pairPotential().scalarPotentialForce(distance, interactionType);

    particle1.force() -= force12 * r12;
    particle2.force() += force12 * r12;
  }
  
  ///Computes the force and the energy associated to this pair interaction, and updates these 2 fields
  ///The first 2 derivates of the potential are stored in "particle2"
  ///Here the distance is not periodic !
  void Colloid::nonPerInteraction(Particle& particle1, Particle& particle2) const
  {
    int interactionType = particle1.type() + particle2.type();
    
    DVec r12 = particle2.position() - particle1.position();
    double distance = r12.norm();
    r12 /= distance;
    // compute energy
    double energy12 = pairPotential()(distance, interactionType);

    particle1.potentialEnergy() += energy12 / 2;
    particle2.potentialEnergy() += energy12 / 2;
    // compute forces (positive -> repulsive / negative -> attractive)
    double force12 = pairPotential().scalarPotentialForce(distance, interactionType);
    particle1.force() -= force12 * r12;
    particle2.force() += force12 * r12;
  }
  
  void Colloid::samplePositions(DynamicsParameters const& dynaPara)
  {
    NBody::samplePositions(dynaPara);
    
    //!\ Can put the system in an unstable configuration if an end of the dimer starts too close to a particle of solvent
    if (dynaPara.interactionRatio() == 0)
    {
      DVec q0 = getParticle(0).position();
      DVec q1 = getParticle(1).position();
      DVec r10 = q1 - q0;
      r10 /= r10.norm();
      DVec qMid = (q0+q1)/2;
      double dist = pairPotential_->drawLaw(dynaPara.beta(), rng_, 1);
      getParticle(0).position() = qMid - dist/2 * r10;
      getParticle(1).position() = qMid + dist/2 * r10;
    }
  }
  
  
  
  void Colloid::computeAllForces()
  {
    for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
      getParticle(iOfParticle).resetForce(externalPotential());
     
    if (doCells_)
    {
      //-- reinitialize cells before looping on the pair interactions --
      reinitializeCells();
      for (ParticlePairIterator it = pairBegin(); !pairFinished(it); incrementePairIterator(it))
        interaction(it.particle1(), it.particle2());
      
      // we make the two ends of the colloid interact even if they are not in the same box, and we don't periodize their positions
      nonPerInteraction(getParticle(0), getParticle(1));
    }
    else
    {
      //-- no cell method: std double loop --
      for (int i = 0; i < nbOfParticles(); i++)
        for (int j = i + 1; j < nbOfParticles(); j++)
          interaction(getParticle(i), getParticle(j));
        
      nonPerInteraction(getParticle(0), getParticle(1));
    } 
  }
  
  ///
  ///Computes the instant value of the observable length
  ///This length is not periodized!
  double Colloid::length() const
  {
    DVec r01 = getParticle(1).position() - getParticle(0).position();
    return r01.norm();
  }
  
  ///
  ///Computes the instant value of the observable force
  /// positive -> repulsion / negative -> attraction
  double Colloid::force() const
  {
    DVec r01 = getParticle(1).position() - getParticle(0).position();
    double d01 = r01.norm();
   
    return dot(getParticle(1).force() - getParticle(0).force(), r01) / (2*d01);
  }
  
}