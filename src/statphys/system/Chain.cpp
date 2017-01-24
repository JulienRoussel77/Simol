//#include <core/io/CommandLine.hpp>

#include "simol/statphys/system/Chain.hpp"

namespace simol
{
  //### Chain ###

  Chain::Chain(Input const& input, int nbOfWallParticles0):
    System(input),
    nbOfWallParticles_(nbOfWallParticles0)
  {
    configuration_ = vector<Particle*>(nbOfParticles()+nbOfWallParticles_, nullptr);    
    for (int iOfParticle = -nbOfWallParticles_; iOfParticle < input.nbOfParticles(); iOfParticle++)
      configuration_[iOfParticle+nbOfWallParticles_] = new Particle(input.mass(), input.initialPosition(iOfParticle), input.initialMomentum(iOfParticle));
    
    
    for (int iOfParticle = -nbOfWallParticles_; iOfParticle < 0; iOfParticle++)
    {
      getParticle(iOfParticle).position(0) = 0;
      getParticle(iOfParticle).momentum(0) = 0;
    }
  }

  //--------------- particle pair iterators ------------
  
  void Chain::incrementePairIterator(ParticlePairIterator& it)
  {
    //it.iOfParticle1()++;
    //it.iOfParticle2()++;
    it.it1_++;
    it.it2_++;
  }

  bool Chain::pairFinished(ParticlePairIterator const& it) const
  {
    //return it.iOfParticle1() == nbOfParticles() - 1;
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
  


  //###### BiChain ######

  BiChain::BiChain(Input const& input):
    Chain(input, 0)
  {}



  void BiChain::samplePositions(DynamicsParameters const& dynaPara)
  {
    cout << " - Sampling the positions..." << endl;
    double alpha, localTemp, localDist;
    for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
    {
      alpha = iOfParticle / (double) nbOfParticles();
      localTemp = (1 - alpha) * dynaPara.temperatureLeft() + alpha * dynaPara.temperatureRight();
      localDist = drawPotLaw(1 / localTemp);
      double prevPosition = (iOfParticle==0)?0:getParticle(iOfParticle - 1).position(0);

      getParticle(iOfParticle).position(0) = prevPosition + localDist;
    }
  }



  void BiChain::computeAllForces()
  {
    //std::cout << "BiChain::computeAllForces" << std::endl;
    //for (auto && particle : configuration_)
    //  particle.resetForce(potential());
    //for (int i = 0; i < nbOfParticles() - 1; i++)
    // interaction(configuration_[i], configuration_[i + 1]);
    
    //for (ParticleIterator it = begin(); !finished(it); ++it)
    for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
      getParticle(iOfParticle).resetForce(potential());
    
    for (ParticlePairIterator it = pairBegin(); !pairFinished(it); incrementePairIterator(it))
      interaction(it.particle1(), it.particle2());
  }




  //###### TriChain ######


  TriChain::TriChain(Input const& input):
    Chain(input, 2),
    /*getParticle(-2)(input.dimension()),
    ancorParticle2_(input.dimension()),*/
    isOfFixedVolum_(input.isOfFixedVolum())
  {}

  bool const& TriChain::isOfFixedVolum() const
  {
    return isOfFixedVolum_;
  }
  
  void TriChain::samplePositions(DynamicsParameters const& dynaPara)
  {
    cout << " - Sampling the positions..." << endl;
    if (isOfFixedVolum())
    {
      cout << "    - Simulation of fixed volum : q_i = 0..." << endl;
      for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
        getParticle(iOfParticle).position(0) = 0;
    }
    else
    {
      double alpha, localTemp, localBending;
      for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
      {
        alpha = iOfParticle / (double) nbOfParticles();
        localTemp = (1 - alpha) * dynaPara.temperatureLeft() + alpha * dynaPara.temperatureRight();

        localBending = drawPotLaw(1 / localTemp);
        double position1 = getParticle(iOfParticle - 1).position(0);
        double position2 = getParticle(iOfParticle - 2).position(0);
        getParticle(iOfParticle).position(0) = -position2 + 2 * position1 + localBending;
        getParticle(iOfParticle).momentum() = drawMomentum(1 / localTemp, getParticle(iOfParticle).mass());
      }
    }
  }

  ///Computes the force and the energy associated to this triplet interaction, and updates these 2 fields
  ///The first 2 derivates of the potential are stored in "particle2"
  void TriChain::triInteraction(Particle& particle1, Particle& particle2, Particle& particle3) const
  {
    DVec delta = particle3.position() - 2 * particle2.position() + particle1.position();
    //double d12 = r12.norm();
    double energy123 = potential(delta);
    DVec force123 = potentialForce(delta);    // = - v'(r_2)
    double lapla123 = laplacian(delta);

    particle2.potentialEnergy() = energy123;
    particle1.force() += force123;
    particle2.force() -= 2 * force123;
    particle3.force() += force123;
    particle2.energyGrad() = -force123;    // - v'(r_2)
    particle2.energyLapla() = lapla123;    // v''(r_2)
  }

  void TriChain::computeAllForces()
  {
    //std::cout << "TriChain::computeAllForces" << std::endl;
    for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
      getParticle(iOfParticle).resetForce(potential());
    //for (auto&& particle : configuration_)
    //  dyna.computeForce(particle);
    triInteraction(getParticle(-2), getParticle(-1), getParticle(0));
    triInteraction(getParticle(-1), getParticle(0), getParticle(1));
    for (int iOfParticle = 0; iOfParticle < nbOfParticles() - 2; iOfParticle++)
      triInteraction(getParticle(iOfParticle), getParticle(iOfParticle + 1), getParticle(iOfParticle + 2));
    //dyna.bending(configuration_[nbOfParticles() - 2], configuration_[nbOfParticles() - 1]);
    if (isOfFixedVolum_)
      triInteraction(getParticle(nbOfParticles() - 1), getParticle(-2), getParticle(-1));
  }

  double TriChain::boundaryPotEnergy() const
  {return getParticle(-2).potentialEnergy() + getParticle(-1).potentialEnergy();}




}
