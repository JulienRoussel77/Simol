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
    //if ((particle1.type() == 1) || (particle2.type() == 1))
    //cout << particle1.type() << " x " << particle2.type() << " -> " << interactionType << endl;
    // take closest periodic image
    //if (particle1.type() != particle2.type()) return;
    DVec r12 = periodicDistance(particle1.position() - particle2.position());
    //DVec r12 = particle1.position() - particle2.position();
    double distance = r12.norm();
    r12 /= distance;
    // compute energy
    double energy12 = pairPotential()(distance, interactionType);

    particle1.potentialEnergy() += energy12 / 2;
    particle2.potentialEnergy() += energy12 / 2;
    // compute forces
    double force12 = pairPotential().scalarPotentialForce(distance, interactionType);
    /*if (interactionType == 0)
    {
      ofstream test("test", std::ofstream::app);
      test << distance << " " << energy12 << " " << force12 << endl;
    }*/
    particle1.force() += force12 * r12;
    particle2.force() -= force12 * r12;
    // compute pressure, based on the Virial formula for the potential part: P_pot = -\sum_{i < j} r_ij v'(r_ij) / d|Vol|; will divide by d|Vol| at the end
    particle1.virial() += 0.5 * force12 * distance;
    particle2.virial() += 0.5 * force12 * distance;
    
    //cout << "inter : " << distance << " " << energy12 << " " << force12 << endl;
    //cout << "--> " << particle1.position().adjoint() << " v " << particle2.position().adjoint() << " -> distance = " << distance << endl;
  }
  
  void Colloid::samplePositions(DynamicsParameters const& dynaPara)
  {
    NBody::samplePositions(dynaPara);
    
    //!\ Can put the system in an unstable configuration...
    if (dynaPara.interactionRatio() == 0)
    {
      DVec q0 = getParticle(0).position();
      DVec q1 = getParticle(1).position();
      DVec r10 = q1 - q0;
      r10 /= r10.norm();
      DVec qMid = (q0+q1)/2;
      double dist = pairPotential_->drawLaw(dynaPara.beta(), rng_, 1);
      cout << "  -> colloid length : " << dist << endl;
      getParticle(0).position() = qMid - dist/2 * r10;
      getParticle(1).position() = qMid + dist/2 * r10;
    }
  }
  
  void Colloid::computeAllForces()
  {
    //cout << endl << "------------------computeAllForces----------------------" << endl;
    //for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
    //  cout << "1 : i " << iOfParticle << " " << getParticle(iOfParticle).type() << endl;

    for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
      getParticle(iOfParticle).resetForce(externalPotential());
    
    //for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
    //  cout << "2 : i " << iOfParticle << ", p " << getParticle(iOfParticle).position().adjoint() << ", e " << getParticle(iOfParticle).potentialEnergy() << ", f " << getParticle(iOfParticle).force().adjoint() << endl;
    
    if (doCells_)
    {
      //-- reinitialize cells before looping on the pair interactions --
      reinitializeCells();
      for (ParticlePairIterator it = pairBegin(); !pairFinished(it); incrementePairIterator(it))
        interaction(it.particle1(), it.particle2());
    }
    else
    {
      //-- no cell method: std double loop --
      for (int i = 0; i < nbOfParticles(); i++)
        for (int j = i + 1; j < nbOfParticles(); j++)
          interaction(getParticle(i), getParticle(j));
    } 

  }
  
  ///
  ///Computes the instant value of the observable length
  double Colloid::length() const
  {
    DVec r12 = periodicDistance(getParticle(0).position() - getParticle(1).position());
    //cout << periodicImage(getParticle(0).position()).adjoint() << " <-> " << periodicImage(getParticle(1).position()).adjoint() << " = " << r12.norm() << endl;
    return r12.norm();
  }
  
}