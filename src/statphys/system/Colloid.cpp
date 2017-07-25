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
    //DVec r12 = particle2.position() - particle1.position();
    double distance = r12.norm();
    r12 /= distance;
    // compute energy
    double energy12 = pairPotential()(distance, interactionType);

    particle1.potentialEnergy() += energy12 / 2;
    particle2.potentialEnergy() += energy12 / 2;
    // compute forces (positive -> repulsive / negative -> attractive)
    double force12 = pairPotential().scalarPotentialForce(distance, interactionType);
    /*if (interactionType == 2)
    {
      ofstream test("test", std::ofstream::app);
      test << distance << " " << energy12 << " " << force12 << endl;
    }*/
    particle1.force() -= force12 * r12;
    particle2.force() += force12 * r12;
    
    //cout << "inter : " << distance << " " << energy12 << " " << force12 << endl;
    //cout << "--> " << particle1.position().adjoint() << " v " << particle2.position().adjoint() << " -> distance = " << distance << endl;
    
    //cout << particle1.type() << " x " << particle2.type() << " : " << distance << " " << force12 << endl;
  }
  
  ///Computes the force and the energy associated to this pair interaction, and updates these 2 fields
  ///The first 2 derivates of the potential are stored in "particle2"
  ///Here the distance is not periodic !
  void Colloid::nonPerInteraction(Particle& particle1, Particle& particle2) const
  {
    int interactionType = particle1.type() + particle2.type();
    
    //DVec r12 = periodicDistance(particle2.position() - particle1.position());
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
      
      nonPerInteraction(getParticle(0), getParticle(1));
      
      // When the two dimer particles are not in neighboring cells we make them interact
      /*int iOfCell0 = findIndex(getParticle(0).position());
      vector<int>* indexNeigh0 = &cell(iOfCell0).indexNeighbors();
      int iOfCell1 = findIndex(getParticle(1).position());
      vector<int>* indexNeigh1 = &cell(iOfCell1).indexNeighbors();
      if (iOfCell0 != iOfCell1 
        && std::find(indexNeigh0->begin(), indexNeigh0->end(), iOfCell1) == indexNeigh0->end()
        && std::find(indexNeigh1->begin(), indexNeigh1->end(), iOfCell0) == indexNeigh1->end())
        interaction(getParticle(0), getParticle(1));*/
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
    //DVec r01 = periodicDistance(getParticle(1).position() - getParticle(0).position());
    DVec r01 = getParticle(1).position() - getParticle(0).position();
    //cout << periodicImage(getParticle(0).position()).adjoint() << " <-> " << periodicImage(getParticle(1).position()).adjoint() << " = " << r12.norm() << endl;
    return r01.norm();
  }
  
  ///
  ///Computes the instant value of the observable force
  /// positive -> repulsion / negative -> attraction
  double Colloid::force() const
  {
    //DVec r01 = periodicDistance(getParticle(1).position() - getParticle(0).position());
    DVec r01 = getParticle(1).position() - getParticle(0).position();
    double d01 = r01.norm();
    //cout << "Out force : " << d01 << " " << dot(getParticle(1).force() - getParticle(0).force(), r01) / (2*d01) << endl;
    return dot(getParticle(1).force() - getParticle(0).force(), r01) / (2*d01);
  }
  
}