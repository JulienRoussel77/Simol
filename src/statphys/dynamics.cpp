#ifndef SIMOL_DYNAMICS_IPP
#define SIMOL_DYNAMICS_IPP

#include "dynamics.hpp"

namespace simol
{

    Dynamics* createDynamics(Input const& input, size_t iOfReplica)
  {
    if (input.dynamicsName() == "Hamiltonian")
      return new Hamiltonian(input, iOfReplica);
    else if (input.dynamicsName() == "Langevin")
      return new Langevin(input, iOfReplica);
    else if (input.dynamicsName() == "BoundaryLangevin")
      return new BoundaryLangevin(input, iOfReplica);
    else if (input.dynamicsName() == "Overdamped")
      return new Overdamped(input, iOfReplica);
    else
      std::cout << input.dynamicsName() << " is not a valid dynamics !" << std::endl;
    return 0;
  }

  Dynamics::Dynamics(Input const& input, int const& iOfReplica):
    timeStep_(input.timeStep(iOfReplica)),
    numberOfIterations_(input.numberOfIterations(iOfReplica)),
    externalForce_(input.dimension())
  {
    potential_ = createPotential(input);
    externalForce_(0) = input.externalForce(iOfReplica);
    std::cout << "externalForce = " << externalForce_(0) << std::endl;
    std::cout << "timeStep = " << timeStep() << std::endl;
    if (numberOfIterations_ < 1e6)
      std::cout << "numberOfIterations_ = " << numberOfIterations_ << std::endl;
    else
      std::cout << "numberOfIterations_ = " << numberOfIterations_/1e6 << "e6" << std::endl;
    std::cout << "finalTime = " << numberOfIterations_ * timeStep() << std::endl;
  }

  Dynamics::~Dynamics()
  {
    delete potential_;
  }

  double& Dynamics::timeStep()
  {
    return timeStep_;
  }

  const double& Dynamics::timeStep() const
  {
    return timeStep_;
  }

  size_t& Dynamics::numberOfIterations()
  {
    return numberOfIterations_;
  }

  const size_t& Dynamics::numberOfIterations() const
  {
    return numberOfIterations_;
  }

  double Dynamics::finalTime() const
  {
    return timeStep_ * numberOfIterations_;
  }

  Potential& Dynamics::potential()
  { return *potential_; }

  double Dynamics::potential(Vector<double> const& position) const
  {
    return (*potential_)(position);
  }

/*  double Dynamics::potential(double const& distance) const
  {
    return (*potential_)(distance);
  }
  */
  Vector<double> Dynamics::force(Vector<double> const& position) const
  {
    return potential_->force(position) + externalForce_;
  }

  const Vector<double>& Dynamics::externalForce() const
  {
    return externalForce_;
  }

  Vector<double>& Dynamics::externalForce()
  {
    return externalForce_;
  }

  const double& Dynamics::externalForce(int const& i) const
  {
    return externalForce_(i);
  }

  double& Dynamics::externalForce(int const& i)
  {
    return externalForce_(i);
  }




  void Dynamics::initializeMomenta(std::vector<Particle>& /*configuration*/)
  {}

  void Dynamics::resetForce(Particle& particle) const
  {
    particle.potentialEnergy() = 0;
    particle.force().fill(0);
  }

  void Dynamics::computeForce(Particle& particle) const
  {
    particle.potentialEnergy() = potential(particle.position());
    particle.force() = force(particle.position());
  }

  void Dynamics::interaction(Particle& particle1, Particle& particle2) const
  {
    Vector<double> r12 = particle2.position() - particle1.position();
    //double d12 = r12.norm();
    double energy12 = potential(r12);
    Vector<double> force12 = force(r12);    // = - v'(q_2 - q_1)

    //particle1.potentialEnergy() += energy12 / 2;
    particle2.potentialEnergy() = energy12;
    particle1.force() -= force12;
    particle2.force() += force12;
    particle1.energyGrad() = -force12;    // v'(q_2 - q_1)
  }

    void Dynamics::triInteraction(Particle& particle1, Particle& particle2, Particle& particle3) const
  {
    Vector<double> delta = particle3.position() - 2*particle2.position() + particle1.position();
    //double d12 = r12.norm();
    double energy123 = potential(delta);
    Vector<double> force123 = force(delta);    // = - v'(r_2)

    //particle1.potentialEnergy() += energy12 / 2;
    particle2.potentialEnergy() = energy123;
    particle1.force() += force123;
    particle2.force() -= 2*force123;
    particle3.force() += force123;
    particle3.energyGrad() = -force123;    // - v'(r_2)
  }

  void Dynamics::updateBefore(Particle& particle)
  {
    //std::cout << "Dynamics::updateBefore" <<std::endl;
    particle.momentum() += timeStep_ * particle.force() / 2;
    particle.position() += timeStep_ * particle.momentum() / particle.mass();
  }

  void Dynamics::updateAfter(Particle& particle)
  {
    //std::cout << "Dynamics::updateAfter" <<std::endl;
    particle.momentum() += timeStep_ * particle.force() / 2;
    particle.kineticEnergy() = pow(particle.momentum().norm(), 2) / particle.mass() / 2;
  }

  void Dynamics::updateAllControlVariates(Output& output, std::vector<Particle> const& configuration, size_t indexOfIteration) const
  {
    Vector<double> q = configuration[0].position();
    Vector<double> p = configuration[0].momentum();
    Vector<double> qEnd = configuration[configuration.size()-1].position();
    //assert(configuration.size() == 1);
    VectorXd generatorOnBasis;
    generatorOnBasis = generatorOn(output.velocityCV(), configuration);
    output.velocityCV()->update(p(0), generatorOnBasis, configuration, indexOfIteration);
    generatorOnBasis = generatorOn(output.forceCV(), configuration);
    output.forceCV()->update(potential_->derivative(q)(0), generatorOnBasis, configuration, indexOfIteration);
    generatorOnBasis = generatorOn(output.lengthCV(), configuration);
    output.lengthCV()->update(qEnd(0), generatorOnBasis, configuration, indexOfIteration);
    generatorOnBasis = generatorOn(output.midFlowCV(), configuration);
    output.midFlowCV()->update(output.energyMidFlow(), generatorOnBasis, configuration, indexOfIteration);
		generatorOnBasis = generatorOn(output.sumFlowCV(), configuration);
    output.sumFlowCV()->update(output.energySumFlow(), generatorOnBasis, configuration, indexOfIteration);
  }

  double Dynamics::computeMeanPotLaw(double localBeta) const
	{
		//std::cout << "Dynamics::computeMeanPotLaw" << std::endl;
		double repFunc = 0;
		double qInteg = 0;
		size_t nbIntegrationNodes = 1000;
		double step = 8. / nbIntegrationNodes;
		Vector<double> deltaQ(1);
		for (size_t iOfNode = 0; iOfNode < nbIntegrationNodes; iOfNode++)
		{
			deltaQ(0) = - 4 + iOfNode * step;
			//std::cout << q << " " << localBeta << " " << potential(q) << std::endl;
			repFunc += exp(-localBeta * potential(deltaQ));
			qInteg += deltaQ(0) * exp(-localBeta * potential(deltaQ));
		}
		//std::cout << "beta = " << localBeta << " -> " << qInteg / repFunc << std::endl;
		return qInteg / repFunc;
	}


  //#### Hamiltonian ####

  Hamiltonian::Hamiltonian(Input const& input, int const& iOfReplica):Dynamics(input, iOfReplica)
  {
  }

  MatrixXd Hamiltonian::generatorOn(ControlVariate const* controlVariate, std::vector<Particle> const& configuration) const
  {
    VectorXd result = VectorXd::Zero(controlVariate->nbOfFunctions());
    for (size_t iOfFunction=0; iOfFunction < controlVariate->nbOfFunctions(); iOfFunction++)
      for (size_t iOfParticle=0; iOfParticle < configuration.size(); iOfParticle++)
				result(iOfFunction) += configuration[iOfParticle].momentum().dot(controlVariate->gradientQ(configuration, iOfParticle, iOfFunction))
	      + force(configuration[iOfParticle].position()).dot(controlVariate->gradientP(configuration, iOfParticle, iOfFunction));
    //return momentum.dot(controlVariate->gradientQ(configuration, i))
    //+ force(position).dot(controlVariate->gradientP(configuration));
    return result;
  }



    //#### StochasticDynamics ####

  StochasticDynamics::StochasticDynamics(Input const& input, int const& iOfReplica):
    Dynamics(input, iOfReplica)
  {}

  void StochasticDynamics::setRNG(RNG* rng)
  {
    rng_ = rng;
  }

  Vector<double> StochasticDynamics::drawMomentum(double localBeta, double mass)
	{
		return sqrt(1 / (localBeta * mass)) * rng_->gaussian();
	}



  //#### UniformStochasticDynamics ####

  UniformStochasticDynamics::UniformStochasticDynamics(Input const& input, int const& iOfReplica):
    StochasticDynamics(input, iOfReplica),
    beta_(input.beta(iOfReplica)),
    temperature_(1/beta_)
  {}


  void UniformStochasticDynamics::initializeMomenta(std::vector<Particle>& configuration)
  {
    for (auto&& particle : configuration)
      particle.momentum() = sqrt(1 / (beta() * particle.mass())) * rng_->gaussian();
  }

  double UniformStochasticDynamics::temperature() const
  {
    return temperature_;
  }

  double const& UniformStochasticDynamics::beta() const
  {
    return beta_;
  }

  double const& UniformStochasticDynamics::temperatureLeft() const
  {
    return temperature_;
  }

  double const& UniformStochasticDynamics::temperatureRight() const
  {
    return temperature_;
  }

  //#### Langevin ####

  Langevin::Langevin(Input const& input, int const& iOfReplica):
    UniformStochasticDynamics(input, iOfReplica),
    gamma_(input.gamma())
  {}


  double const& Langevin::gamma() const
  {
    return gamma_;
  }

  double Langevin::sigma() const
  {
    return sqrt(2*gamma_ / beta_);
  }


  void Langevin::updateAfter(Particle& particle)
  {
    particle.momentum() += timeStep_ * particle.force() / 2;
    double alpha = exp(- gamma_ / particle.mass() * timeStep_);
    particle.momentum() = alpha * particle.momentum() + sqrt((1-pow(alpha, 2))/beta_*particle.mass()) * rng_->gaussian();

    particle.kineticEnergy() = pow(particle.momentum().norm(), 2) / particle.mass() / 2;
  }

  MatrixXd Langevin::generatorOn(ControlVariate const* controlVariate, std::vector<Particle> const& configuration) const
  {
    VectorXd result = VectorXd::Zero(controlVariate->nbOfFunctions());
    for (size_t iOfFunction=0; iOfFunction < controlVariate->nbOfFunctions(); iOfFunction++)
      for (size_t iOfParticle=0; iOfParticle < configuration.size(); iOfParticle++)
	result(iOfFunction) += dot(configuration[iOfParticle].momentum(), controlVariate->gradientQ(configuration, iOfParticle, iOfFunction))
	      + dot(configuration[iOfParticle].force(), controlVariate->gradientP(configuration, iOfParticle, iOfFunction))
	      + gamma_ * (- dot(configuration[iOfParticle].momentum(), controlVariate->gradientP(configuration, iOfParticle, iOfFunction))
			+ controlVariate->laplacienP(configuration, iOfParticle, iOfFunction) / beta_ );
    return result;
    //return configuration[0].position(0);
  }




  //#### Overdamped ####

  Overdamped::Overdamped(Input const& input, int const& iOfReplica):
    UniformStochasticDynamics(input, iOfReplica)
  {}


  void Overdamped::updateBefore(Particle& /*particle*/)
  {

  }

  // /!\ a changer dans le cas N != 1
  void Overdamped::updateAfter(Particle& particle)
  {
    particle.position() += timeStep_ * particle.force() + sqrt(2*timeStep_/beta_) * rng_->gaussian();
  }

  MatrixXd Overdamped::generatorOn(ControlVariate const* controlVariate, std::vector<Particle> const& configuration) const
  {
    VectorXd result = VectorXd::Zero(controlVariate->nbOfFunctions());
      for (size_t iOfFunction=0; iOfFunction < controlVariate->nbOfFunctions(); iOfFunction++)
	for (size_t iOfParticle=0; iOfParticle < configuration.size(); iOfParticle++)
	  result(iOfFunction) += controlVariate->laplacienQ(configuration, iOfParticle, iOfFunction) / beta_
	      + dot(force(configuration[iOfParticle].position()), controlVariate->gradientQ(configuration, iOfParticle, iOfFunction));
    return result;

    //double q = configuration[0].position(0);
    //return - pow (2 * M_PI, 2) * sin(2 * M_PI * q) / beta_ + force(configuration[0].position())(0) * 2 * M_PI * cos(2 * M_PI * q);
  }


  //#### BoundaryStochasticDynamics ####

  BoundaryStochasticDynamics::BoundaryStochasticDynamics(Input const& input, int const& iOfReplica):
    StochasticDynamics(input, iOfReplica),
    betaLeft_(input.betaLeft()),
    betaRight_(input.betaRight()),
    temperatureLeft_(1/betaLeft_),
    temperatureRight_(1/betaRight())
  {
      std::cout << "deltaTemperature = " << deltaTemperature() << std::endl;
  }

    double const& BoundaryStochasticDynamics::betaLeft() const
  {
    return betaLeft_;
  }

  double const& BoundaryStochasticDynamics::betaRight() const
  {
    return betaRight_;
  }

  double BoundaryStochasticDynamics::temperature() const
  {
    return (temperatureLeft_ + temperatureRight_) / 2;
  }

  double const& BoundaryStochasticDynamics::temperatureLeft() const
  {
    return temperatureLeft_;
  }

  double const& BoundaryStochasticDynamics::temperatureRight() const
  {
    return temperatureRight_;
  }

  double BoundaryStochasticDynamics::deltaTemperature() const
  {
    return (temperatureLeft_ - temperatureRight_)/2;
  }

  void BoundaryStochasticDynamics::initializeMomenta(std::vector<Particle>& configuration)
  {
    for (size_t iOfParticle = 0; iOfParticle < configuration.size(); iOfParticle++)
    {
      Particle& particle = configuration[iOfParticle];
      double tempi_ = temperatureLeft_ + (iOfParticle + .5) * (temperatureRight_ - temperatureLeft_) / configuration.size();
      particle.momentum() = sqrt(tempi_ / particle.mass()) * rng_->gaussian();
    }
  }

  /*void BoundaryStochasticDynamics::startFrom(string settingsPath)
  {
    for (size_t iOfParticle = 0; iOfParticle < configuration.size(); iOfParticle++)
    {
      Particle& particle = configuration[iOfParticle];
      double tempi_ = temperatureLeft_ + (iOfParticle + .5) * (temperatureRight_ - temperatureLeft_) / configuration.size();
      particle.momentum() = sqrt(tempi_ / particle.mass()) * rng_->gaussian();
    }
  }*/


  //#### BoundaryLangevin ####

  BoundaryLangevin::BoundaryLangevin(const Input& input, const int& iOfReplica):
    BoundaryStochasticDynamics(input, iOfReplica),
    gamma_(input.gamma()),
    tauBending_(input.tauBending())
  {}

  double const& BoundaryLangevin::gamma() const
  {
    return gamma_;
  }

  double BoundaryLangevin::sigmaLeft() const
  {
    return sqrt(2 * gamma_ / betaLeft_);
  }

  double BoundaryLangevin::sigmaRight() const
  {
    return sqrt(2 * gamma_ / betaRight_);
  }

  double const& BoundaryLangevin::tauBending() const
  {
    return tauBending_;
  }

    void BoundaryLangevin::updateAfterLeft(Particle& particle)
  {
    //particle.momentum() += timeStep_ * particle.force() / 2;
    double alpha = exp(- gamma_ / particle.mass() * timeStep_);
    particle.momentum() = alpha * particle.momentum() + sqrt((1-pow(alpha, 2))/betaLeft_*particle.mass()) * rng_->gaussian();

    particle.kineticEnergy() = pow(particle.momentum().norm(), 2) / particle.mass() / 2;
  }

    void BoundaryLangevin::updateAfterRight(Particle& particle)
  {
    //particle.momentum() += timeStep_ * particle.force() / 2;
    double alpha = exp(- gamma_ / particle.mass() * timeStep_);
    particle.momentum() = alpha * particle.momentum() + sqrt((1-pow(alpha, 2))/betaRight_*particle.mass()) * rng_->gaussian();

    particle.kineticEnergy() = pow(particle.momentum().norm(), 2) / particle.mass() / 2;
  }

  void BoundaryLangevin::bending(Particle& particle1, Particle& particle2) const
  {
    particle1.force(0) -= tauBending();
    particle2.force(0) += tauBending();
  }

  MatrixXd BoundaryLangevin::generatorOn(ControlVariate const* controlVariate, std::vector<Particle> const& configuration) const
  {
    size_t numberOfParticles = configuration.size();
    VectorXd result = VectorXd::Zero(controlVariate->nbOfFunctions());
    for (size_t iOfFunction=0; iOfFunction < controlVariate->nbOfFunctions(); iOfFunction++)
    {
      for (size_t iOfParticle=0; iOfParticle < configuration.size(); iOfParticle++)
      result(iOfFunction) += dot(configuration[iOfParticle].momentum(), controlVariate->gradientQ(configuration, iOfParticle, iOfFunction))
	      + dot(configuration[iOfParticle].force(), controlVariate->gradientP(configuration, iOfParticle, iOfFunction));
	      //if(false)
      result(iOfFunction) += gamma_ * (- dot(configuration[0].momentum(), controlVariate->gradientP(configuration, 0, iOfFunction))
			+ controlVariate->laplacienP(configuration, 0, iOfFunction) / betaLeft_
			- dot(configuration[numberOfParticles-1].momentum(), controlVariate->gradientP(configuration, numberOfParticles-1, iOfFunction))
			+ controlVariate->laplacienP(configuration, numberOfParticles-1, iOfFunction) / betaRight_);
    }
    return result;
  }


}

#endif
