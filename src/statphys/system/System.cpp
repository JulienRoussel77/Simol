#ifndef SIMOL_PARTICLESYSTEM_IMPL_HPP
#define SIMOL_PARTICLESYSTEM_IMPL_HPP

#include "simol/statphys/system/System.hpp"

namespace simol
{ 

  
  ///
  /// Assumes that there are 2 particles in the first box ! Can be improved...
  ParticlePairIterator::ParticlePairIterator():
    //syst_(&syst0),
    //it1_(syst0.cell(0).members().begin()),
    iOfCell1_(0),
    //it2_(std::next(it1_, 1)),
    iOfNeighbor2_(0)
  {} 
  
  System::System(Input const& input):
    dimension_(input.dimension()),
    nbOfParticles_(input.nbOfParticles()),
    //configuration_(input.nbOfParticles(), new Particle(dimension_)),
    configuration_(nbOfParticles_, nullptr),
    doSetting_(input.doSetting()),
    domainSize_(std::numeric_limits<double>::infinity()),
    systemSubtype_(input.systemSubtype())
  {
    externalPotential_ = createPotential(input, input.externalPotentialName());
    pairPotential_ = createPotential(input, input.pairPotentialName());
    if (externalPotential_) domainSize_ = externalPotential_->domainSize();
    
    cout << "System::System" << endl;
    cout << "externalPotential_ : " << externalPotential_->classname() << endl;
    cout << "pairPotential_ : " << pairPotential_->classname() << endl;
  }

  ///
  ///Destructor
  System::~System()
  {
    if (externalPotential_) delete externalPotential_;
    if (pairPotential_) delete pairPotential_;
    for (int iOfParticle = 0; iOfParticle < (int)configuration_.size(); iOfParticle++)
      delete configuration_[iOfParticle];
  }

  void System::printName() const
  {
    cout << "SystemType = System" << endl;
  }

  const Particle& System::getParticle(int index) const
  {
    return operator()(index);
  }

  Particle& System::getParticle(int index)
  {
    return operator()(index);
  }

  const int& System::dimension() const
  {
    return dimension_;
  }

  const std::vector<Particle*> & System::configuration() const
  {
    return configuration_;
  }

  std::vector<Particle*> & System::configuration()
  { return configuration_; }

  const int& System::nbOfParticles() const
  {
    return nbOfParticles_;
  }
  
  const string& System::systemSubtype() const
  {
    return systemSubtype_;
  }
  
  ///
  /// return position x in the interval [0,domainSize_)
  double System::periodicImage(double x) const
  {
    return x - floor(x / domainSize_) * domainSize_;
  }
  
  ///
  /// return position x in the domain [0,domainSize_)^dim
  DVec System::periodicImage(DVec const& position) const
  {
    if (domainSize_ == std::numeric_limits<double>::infinity()) return position;
    DVec periodicPos = position;
    for (int iOfDim = 0; iOfDim < dimension(); iOfDim++)
      periodicPos[iOfDim] = periodicImage(position[iOfDim]);
    return periodicPos;
  }
  
  DVec System::periodicDistance(DVec const& vecDistance) const
  {
    return vecDistance - rint(vecDistance / domainSize_) * domainSize_;
  }

  //-------------- Random numbers ----------------

  const std::shared_ptr<RNG>& System::rng() const {return rng_;}

  std::shared_ptr<RNG>& System::rng() {return rng_;}
  
  //--------------- particle pair iterators ------------
    
  ParticlePairIterator System::pairBegin()
  {
    ParticlePairIterator iter = ParticlePairIterator();
    iter.it1_ = configuration().begin();
    iter.it2_ = std::next(iter.it1_, 1);
    iter.endIt2_ = configuration().end();
    return iter;
  }
  
  bool System::pairFinished(ParticlePairIterator const& it) const
  {
    return it.it1_ == configuration().end();
  }
  
  void System::incrementePairIterator(ParticlePairIterator& it)
  {
    it.it2_++;
    if (it.it2_ == it.endIt2())
    {
      it.it1_++;
      it.it2_ = std::next(it.it1_, 1);
    }
  }

  //----------------- Potential and forces ---------------
  ///
  ///Returns by value the potential of the dynamics
  Potential& System::externalPotential() {return *externalPotential_;}
  ///
  ///Returns by const value the potential of the dynamics
  Potential const& System::externalPotential() const {return *externalPotential_;}
  ///
  ///Returns by value the potential of the dynamics
  Potential& System::pairPotential() {return *pairPotential_;}
  ///
  ///Returns by const value the potential of the dynamics
  Potential const& System::pairPotential() const {return *pairPotential_;}
  ///
  ///Evaluate the potential for the vector "position"
  double System::externalPotential(DVec const& position, int type) const {return (*externalPotential_)(position, type);}
  ///
  ///Evaluate the potential for the scalar "position"
  double System::externalPotential(const double& distance, int type) const {return (*externalPotential_)(distance, type);}
  /*///
  ///Evaluate the force for the scalar "position" (potential and external terms)
  DVec System::totalForce(DVec const& position, int type) const
  {
    return potential_->totalForce(position, type);
  }*/

  /*DVec System::potentialForce(DVec const& position, int type) const
  {
    return - potential_->gradient(position, type);
  }

  DVec System::potentialForce(double position, int type) const
  {
    return - potential_->gradient(position, type);
  }*/
  
  /*///
  ///Evaluate the force for the scalar "position" (potential and external terms)
  DVec System::totalForce(double position, int type) const
  {
    return potential_->totalForce(position, type);
  }*/
  
  /*///
  ///Evaluate the laplacian of the potential for the vector "position"
  double System::laplacian(DVec const& position, int type) const
  {
    return potential_->laplacian(position, type);
  }*/

  /*///
  ///Read-only accessor for the external force
  DVec const& System::externalForce() const {return externalPotential_->force();}
  ///
  ///Write-read accessor for the external force
  DVec& System::externalForce() {return externalPotential_->force();}*/
  /*///
  ///Read-only accessor for the i-th component of the external force
  double const& System::externalForce(int const& i) const {return potential_->externalForce(i);}
  ///
  ///Write-read accessor for the i-th component of the external force
  double& System::externalForce(int const& i) {return potential_->externalForce(i);}*/
  
  ///
  ///Read-only accessor for the first parameter of the potential
  double const& System::potParameter1() const {return pairPotential_->parameter1();}
  
  ///
  ///Read-only accessor for the first parameter of the potential
  double const& System::potParameter2() const {return pairPotential_->parameter2();}
  
  ///
  /// Returns a matrix containing the force vectors of all the particles
  DMat System::forces() const
  {
    DMat F = DMat::Zero(dimension(), nbOfParticles());
    for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
      F.col(iOfParticle) = getParticle(iOfParticle).force();
    return F;
  }
  
  ///
  /// Returns a matrix containing the momenta vectors of all the particles
  DMat System::momenta() const
  {
    DMat P = DMat::Zero(dimension(), nbOfParticles());
    for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
      P.col(iOfParticle) = getParticle(iOfParticle).momentum();
    return P;
  }

  ///
  ///Draw a momentum under the invariant measure at inverse temperature "localBeta"
  DVec System::drawMomentum(double localBeta, double mass)
  {
    return sqrt(1 / (localBeta * mass)) * rng_->gaussian();
  }
  ///
  ///Draw a distance or a bending under the invariant measure at inverse temperature "localBeta"
  double System::drawPotLaw(double localBeta)
  {
    return externalPotential_->drawLaw(localBeta, rng_);
  }
  
  void System::samplePositions(DynamicsParameters const&)
  {
    printName();
    throw std::invalid_argument("samplePositions : Function undefined");
  }
  
  //-- initialization of the momenta according to a Gaussian distribution --
  void System::sampleMomenta(DynamicsParameters const& dynaPara)
  {
    cout << " - Sampling the momenta..." << endl;
    for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
      getParticle(iOfParticle).momentum() = drawMomentum(dynaPara.beta(), getParticle(iOfParticle).mass());
  }

  void System::computeAllForces()
  {throw std::invalid_argument("compteAllForces not defined");}
  

  //--------------- move into chain -----------------

  double System::boundaryPotEnergy() const
  {return 0;}


  ///Computes the force and the energy associated to this pair interaction, and updates these 2 fields
  ///The first 2 derivates of the potential are stored in "particle2"
  void System::interaction(Particle& particle1, Particle& particle2) const
  {
    DVec r12 = particle2.position() - particle1.position();
    double energy12 = pairPotential()(r12);
    DVec force12 = pairPotential_->potentialForce(r12);    // = - v'(q_2 - q_1)
    double lapla12 = pairPotential_->laplacian(r12);  // v"(q_2 - q_1)

    particle1.potentialEnergy() = energy12;
    particle1.force() -= force12;
    particle2.force() += force12;
    particle1.energyGrad() = -force12;    // v'(q_2 - q_1)
    particle1.energyLapla() = lapla12;    // v"(q_2 - q_1)
  }

  ///Compute the mean distance or bending under the invariant measure
  ///Proceeds to a simple integral quadrature using rectangles
  double System::computeMeanPotLaw(double localBeta) const
  {
    double repFunc = 0;
    double qInteg = 0;
    int nbIntegrationNodes = 1000;
    double step = 8. / nbIntegrationNodes;
    DVec deltaQ(1);
    for (int iOfNode = 0; iOfNode < nbIntegrationNodes; iOfNode++)
    {
      deltaQ(0) = - 4 + iOfNode * step;
      repFunc += exp(-localBeta * externalPotential(deltaQ));
      qInteg += deltaQ(0) * exp(-localBeta * externalPotential(deltaQ));
    }
    return qInteg / repFunc;
  }
  
  ///
  ///Computes the instant value of the observable length
  double System::length() const
  {
    return getParticle(0).position(0);
  }
  
  ///
  ///Computes the instant value of the observable velocity
  double System::velocity() const
  {
    return getParticle(0).velocity(0);
  }

  ///
  ///Computes the instant value of the observable force
  double System::force() const
  {
    return getParticle(0).force(0);
  }
}

#endif
