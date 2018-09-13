//#include <core/io/CommandLine.hpp>

#include "simol/system/Chain.hpp"

namespace simol
{
  //### Chain ###

  Chain::Chain(Input const& input):
    System(input)
  {
    configuration_ = vector<Particle*>(nbOfParticles(), nullptr);    
    for (int iOfParticle = 0; iOfParticle < input.nbOfParticles(); iOfParticle++)
      configuration_[iOfParticle] = new Particle(input.mass(), input.initialPosition(iOfParticle), input.initialMomentum(iOfParticle));
    
    
    for (int iOfParticle = 0; iOfParticle < 0; iOfParticle++)
    {
      getParticle(iOfParticle).position(0) = 0;
      getParticle(iOfParticle).momentum(0) = 0;
    }
  }

  //--------------- particle pair iterators ------------
  
  void Chain::incrementePairIterator(ParticlePairIterator& it)
  {
    it.it1_++;
    it.it2_++;
  }

  bool Chain::pairFinished(ParticlePairIterator const& it) const
  {
    return it.it2_ == configuration().end();
  }
  
  void Chain::sampleMomenta(DynamicsParameters const& dynaPara)
  {
    cout << " - Sampling the momenta..." << endl;
    double alpha, localTemp;
    for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
    {
      alpha = iOfParticle / (double) nbOfParticles();
      localTemp = (1 - alpha) * dynaPara.temperatureLeft() + alpha * dynaPara.temperatureRight();
      getParticle(iOfParticle).momentum() = drawMomentum(1 / localTemp, getParticle(iOfParticle).mass());
    }
  }
  
  ///
  ///Computes the instant value of the observable length
  double Chain::length() const
  {
    // /!\ Warning : the chain length is computed assuming it is not attached...
    DVec r12 = getParticle(0).position() - getParticle(nbOfParticles()-1).position();
    return r12.norm();
  }
  
  ///
  ///Draw a distance or a bending under the invariant measure at inverse temperature "localBeta"
  double Chain::drawPotLaw(double localBeta)
  {
    return pairPotential_->drawLaw(localBeta, rng_);
  }
  


  /// The initial positions are sampled by drawing r_i ~ exp(-v(r_i)/T_i)
  /// The temperature T_i at the site i is approximated by a linear profile in the nonequilibrium case
  void Chain::samplePositions(DynamicsParameters const& dynaPara)
  {
    cout << " - Sampling the positions..." << endl;
    double alpha, localTemp, localDist;
    double refPosition = 0;
    for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
    {
      alpha = iOfParticle / (double) nbOfParticles();
      localTemp = (1 - alpha) * dynaPara.temperatureLeft() + alpha * dynaPara.temperatureRight();
      localDist = drawPotLaw(1 / localTemp);

      refPosition += localDist;
      getParticle(iOfParticle).position(0) = refPosition;
      getParticle(iOfParticle).bulkDriving() = dynaPara.bulkDriving();
    }
    getParticle(0).bulkDriving() = 0;
    getParticle(nbOfParticles()-1).bulkDriving() = 0;
  }



  void Chain::computeAllForces()
  {        
    getParticle(0).force(0) = 0; // Only reset the first particle to gain time
    
    for (ParticlePairIterator it = pairBegin(); !pairFinished(it); incrementePairIterator(it))
      interaction(it.particle1(), it.particle2());
  }


  ///Computes the force and the energy associated to this pair interaction, and updates these 2 fields
  ///The first 2 derivates of the potential are stored in "particle2"
  void Chain::interaction(Particle& particle1, Particle& particle2) const
  {
    double r12 = particle2.position(0) - particle1.position(0);
    double energy12 = pairPotential()(r12);
    double force12 = pairPotential_->scalarPotentialForce(r12);    // = - v'(q_2 - q_1)
    double lapla12 = pairPotential_->laplacian(r12);  // v"(q_2 - q_1)

    particle1.potentialEnergy() = energy12;
    particle1.force(0) -= (1-particle1.bulkDriving())*force12;
    particle2.force(0) = (1+particle1.bulkDriving())*force12;    // /!\ Becareful here we assume that the particles are visited in a certain order
    particle1.energyGrad(0) = -force12;    // v'(q_2 - q_1)
    particle1.energyLapla() = lapla12;    // v"(q_2 - q_1)
  }
  
  double Chain::leftHeatFlux(int iOfLink) const
  {
    return - getParticle(iOfLink).momentum(0) / getParticle(iOfLink).mass() * getParticle(iOfLink).energyGrad(0);
  }
  
  double Chain::rightHeatFlux(int iOfLink) const
  {
    return - getParticle(iOfLink+1).momentum(0) / getParticle(iOfLink+1).mass() * getParticle(iOfLink).energyGrad(0);
  }
  
  double Chain::heatFlux(int iOfLink) const
  {
    return (leftHeatFlux(iOfLink) + rightHeatFlux(iOfLink))/2;
  }
  
  double Chain::heatFluxOnAtom(int iOfParticle) const
  {
    return (rightHeatFlux(iOfParticle-1) + leftHeatFlux(iOfParticle))/2;
  }
  
  double Chain::computeSumFlux() const
  {    
    double sumFlux = 0;  
    
    for (int iOfParticle = 1; iOfParticle < nbOfParticles()-1; iOfParticle++)
      sumFlux += heatFluxOnAtom(iOfParticle);      
    return sumFlux / (nbOfParticles() - 2.);
  }
  
  /// Here we compute the sensitivity G(r,p) and project on the right hyperplan
  void Chain::enforceConstraint(double flux, DynamicsParameters const& /*dynaPara*/, bool updateLagrangeMultiplier)
  {    
    double instantFlux = computeSumFlux();
   
    double sensitivity = 0;
    for (int iOfParticle = 1; iOfParticle < nbOfParticles()-1; iOfParticle++)
      sensitivity += pow(getParticle(iOfParticle-1).energyGrad(0) + getParticle(iOfParticle).energyGrad(0), 2);
    
    sensitivity /= pow(2*(nbOfParticles()-2), 2)*getParticle(0).mass();
         
    double localLagrangeMultiplier = -(instantFlux - flux) / sensitivity;

    for (int iOfParticle = 1; iOfParticle < nbOfParticles()-1; iOfParticle++)
      getParticle(iOfParticle).momentum(0) -= localLagrangeMultiplier / (2*(nbOfParticles()-2)) * (getParticle(iOfParticle-1).energyGrad(0) + getParticle(iOfParticle).energyGrad(0));
    
    if (updateLagrangeMultiplier)
      lagrangeMultiplier() += localLagrangeMultiplier;
  }
}
