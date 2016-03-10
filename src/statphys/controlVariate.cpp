#include "controlVariate.hpp"

namespace simol
{

  ControlVariate* createControlVariate(Input const& input, Potential& potential, size_t indexOfReplica)
  {
    if (input.controlVariateName() == "None")
      return new NoControlVariate(input, potential, indexOfReplica);
    if (input.controlVariateName() == "Sinus")
      return new SinusControlVariate(input, potential, indexOfReplica);
    if (input.controlVariateName() == "Cosine")
      return new CosControlVariate(input, potential, indexOfReplica);
    else if (input.controlVariateName() == "SinExp")
      return new SinExpControlVariate(input, potential, indexOfReplica);
    else if (input.controlVariateName() == "CosExp")
      return new CosExpControlVariate(input, potential, indexOfReplica);
    else if (input.controlVariateName() == "Langevin")
      return new LangevinControlVariate(input, potential, indexOfReplica);
    else if (input.controlVariateName() == "SumEnergy")
      return new SumEnergyControlVariate(input, potential, indexOfReplica);
    else if (input.controlVariateName() == "Energy")
      return new EnergyControlVariate(input, potential, indexOfReplica);
    else if (input.controlVariateName() == "Local")
      return new LocalControlVariate(input, potential, indexOfReplica);
    else if (input.controlVariateName() == "Kinetic")
      return new KineticControlVariate(input, potential, indexOfReplica);
    else if (input.controlVariateName() == "Two")
      return new TwoControlVariate(input, potential, indexOfReplica);
    else if (input.controlVariateName() == "IO")
      return new IOControlVariate(input, potential, indexOfReplica);
    else
      std::cout << input.controlVariateName() << " is not a valid control variate !" << std::endl;
    return 0;
  }


  ControlVariate::ControlVariate(Input const& input, Potential& potential, size_t indexOfReplica, size_t nbOfFunctions):
    decorrelationNumberOfIterations_(input.decorrelationNumberOfIterations(indexOfReplica)),
    decorrelationTime_(input.decorrelationTime(indexOfReplica)),
    nbOfFunctions_(nbOfFunctions),
    nbOfFunctionPairs_(pow(nbOfFunctions_, 2)),
    periodNumberOfIterations_(input.outputPeriodNumberOfIterations()),
    statsObservable_(decorrelationNumberOfIterations(), decorrelationTime()),
    statsBetterObservable_(decorrelationNumberOfIterations(), decorrelationTime()),
    statsGeneratorOnBasis_(nbOfFunctions_),
    statsB1_(nbOfFunctions_),
    statsB2_(decorrelationNumberOfIterations(), decorrelationTime(), nbOfFunctions_),
    statsD_(nbOfFunctions_, nbOfFunctions_),
    lastA_(nbOfFunctions_),
    statsPostObservable_(decorrelationNumberOfIterations(), decorrelationTime()),
    statsPostBetterObservable_(decorrelationNumberOfIterations(), decorrelationTime()),
    historyObservable_(input.numberOfIterations(indexOfReplica)),
    historyGeneratorOnBasis_(input.numberOfIterations(indexOfReplica), nbOfFunctions_),
    potential_(&potential)
  {}

  size_t ControlVariate::nbOfFunctions() const
  {
    return nbOfFunctions_;
  }

  size_t ControlVariate::nbOfFunctionPairs() const
  {
    return nbOfFunctionPairs_;
  }

  bool ControlVariate::doOutput(size_t indexOfIteration) const
  {
    return (periodNumberOfIterations_ > 0 && indexOfIteration % periodNumberOfIterations_ == 0);
  }

  bool ControlVariate::isNone() const
  {
    return false;
  }

  double ControlVariate::potential(Vector<double> const& position) const
  {
    return (*potential_)(position);
  }

  Vector<double> ControlVariate::potentialDerivative(Vector<double> const& position) const
  {
    return potential_->derivative(position);
  }

  double ControlVariate::potentialLaplacian(Vector<double> const& position) const
  {
    return potential_->laplacian(position);
  }



  size_t ControlVariate::decorrelationNumberOfIterations() const
  {
    return decorrelationNumberOfIterations_;
  }

  double ControlVariate::decorrelationTime() const
  {
    return decorrelationTime_;
  }



  double ControlVariate::lastObservable() const
  {
    return statsObservable_.lastValue();
  }


  double ControlVariate::meanObservable() const
  {
    return statsObservable_.mean();
  }

  double ControlVariate::stdDeviationObservable() const
  {
    return statsObservable_.standardDeviation();
  }

  double ControlVariate::lastBetterObservable() const
  {
    return statsBetterObservable_.lastValue();
  }

  double ControlVariate::meanBetterObservable() const
  {
    return statsBetterObservable_.mean();
  }

  double ControlVariate::stdDeviationBetterObservable() const
  {
    return statsBetterObservable_.standardDeviation();
  }

  VectorXd ControlVariate::correlationB2() const
  {
    return statsB2_.integratedAutocorrelationMat();
  }

  double ControlVariate::correlationB2(size_t iOfFunction) const
  {
    return statsB2_.integratedAutocorrelation(iOfFunction);
  }

  void ControlVariate::appendToObservable(double observable, size_t iOfIteration)
  {
    historyObservable_(iOfIteration) = observable;
    statsObservable_.append(observable, iOfIteration);
  }

  double ControlVariate::autocorrelation(size_t indexDifference) const
  {
    return statsObservable_(indexDifference);
  }

  void ControlVariate::appendToB1(double observable, VectorXd& valueBasisFunction)
  {
    for (size_t iOfFunction=0; iOfFunction<nbOfFunctions_; iOfFunction++)
      statsB1_.append((observable - statsObservable_.mean()) * valueBasisFunction(iOfFunction), iOfFunction);
  }

    void ControlVariate::appendToB2(double observable, VectorXd& generatorOnBasisFunction, size_t iOfIteration)
  {
    for (size_t iOfFunction=0; iOfFunction<nbOfFunctions_; iOfFunction++)
      statsB2_.append((observable - statsObservable_.mean()), iOfIteration, iOfFunction, generatorOnBasisFunction(iOfFunction));
  }

  void ControlVariate::appendToD(VectorXd& generatorOnBasisFunction, VectorXd& valueBasisFunction)
  {
    //Eigen::Matrix<AutocorrelationStats<double>, Eigen::Dynamic, Eigen::Dynamic> A(nbOfFunctions_, nbOfFunctions_, AutocorrelationStats<double>(decorrelationNumberOfIterations(), decorrelationTime()));

    for (size_t iOfFunction=0; iOfFunction<nbOfFunctions_; iOfFunction++)
      for (size_t iOfFunction2=0; iOfFunction2<=iOfFunction; iOfFunction2++)
      {
	double valueSym = (- valueBasisFunction(iOfFunction) * generatorOnBasisFunction(iOfFunction2)
		      - valueBasisFunction(iOfFunction2) * generatorOnBasisFunction(iOfFunction)) / 2.;
        statsD_.append(valueSym, iOfFunction, iOfFunction2);
	if (iOfFunction != iOfFunction2)
	  statsD_.append(valueSym, iOfFunction2, iOfFunction);
      }
  }



  void ControlVariate::appendToBetterObservable(double observable, VectorXd& generatorOnBasisFunction, size_t iOfIteration)
  {
    //double betterObservableTerm = dot((statsB_.meanMat() - statsB2_.integratedAutocorrelationMat()), statsD_.meanMat().llt().solve(generatorOnBasisFunction));
    /*MatrixXd Dinv = statsD_.meanMat().llt().solve(generatorOnBasisFunction);
    VectorXd B = statsB_.meanMat() - statsB2_.integratedAutocorrelationMat();
    std::cout << B.transpose() * Dinv << std::endl;*/
    //std::cout << (statsB_.meanMat() - statsB2_.integratedAutocorrelationMat()).transpose().size() << "  " << statsD_.meanMat().inverse().size() << std::endl;
    //lastA_ = (statsB_.meanMat() - statsB2_.integratedAutocorrelationMat()).transpose() * statsD_.meanMat().llt().solve(generatorOnBasisFunction);
    //std::cout << lastA_.size() << "   " << generatorOnBasisFunction.size() << std::endl;

    if (statsD_.meanMat().determinant() != 0)
    {
      lastA_ = - .5 * statsD_.meanMat().llt().solve(meanB());
      statsBetterObservable_.append(observable - lastA_.transpose().dot(generatorOnBasisFunction), iOfIteration);
    }
    else
    {
      lastA_.fill(0);
      statsBetterObservable_.append(observable, iOfIteration);
    }

    for (size_t iOfFunction =0; iOfFunction < nbOfFunctions_; iOfFunction++)
    {
      historyGeneratorOnBasis_(iOfIteration, iOfFunction) = generatorOnBasisFunction(iOfFunction);
      statsGeneratorOnBasis_.append(generatorOnBasisFunction(iOfFunction), iOfFunction);
    }
  }



   VectorXd ControlVariate::lastB1() const
  {
    /*VectorXd result(nbOfFunctions_);
    for (size_t iOfFunction =0; iOfFunction < nbOfFunctions_; iOfFunction++)
      result(iOfFunction) = statsB_.lastValue(iOfFunction);
    return result;*/
    return statsB1_.lastValueMat();
  }


  VectorXd ControlVariate::meanB1() const
  {
    return statsB1_.meanMat();
  }

  VectorXd ControlVariate::meanB() const
  {
    return statsB1_.meanMat() - statsB2_.integratedAutocorrelationMat();
  }

  MatrixXd ControlVariate::lastD() const
  {
    return statsD_.lastValueMat();
  }

  MatrixXd ControlVariate::meanD() const
  {
    return statsD_.meanMat();
  }

    VectorXd ControlVariate::lastGeneratorOnBasis() const
  {
    return statsGeneratorOnBasis_.lastValueMat();
  }

  VectorXd ControlVariate::meanGeneratorOnBasis() const
  {
    return statsGeneratorOnBasis_.meanMat();
  }

  double ControlVariate::autocorrelationB2(size_t indexDifference, size_t iOfFunction) const
  {
    return statsB2_(indexDifference, iOfFunction);
  }

  VectorXd ControlVariate::lastA() const
  {
    return lastA_;

  }
  double ControlVariate::lastA(size_t iOfFunction) const
  {
    return lastA_(iOfFunction);
  }


    void ControlVariate::update(double observable, VectorXd& generatorOnBasisFunction, std::vector<Particle> const& configuration, size_t iOfIteration)
  {
    VectorXd valueBasisFunction(nbOfFunctions_);
    for (size_t iOfFunction =0; iOfFunction < nbOfFunctions_; iOfFunction++)
      valueBasisFunction(iOfFunction) = basisFunction(configuration, iOfFunction);
    appendToB1(observable, valueBasisFunction);
    appendToB2(observable, generatorOnBasisFunction, iOfIteration);
    appendToD(generatorOnBasisFunction, valueBasisFunction);
    appendToObservable(observable, iOfIteration);
    appendToBetterObservable(observable, generatorOnBasisFunction, iOfIteration);
  }

    void ControlVariate::display(std::ofstream& out, double time) const
  {
    /*out << time
	<< " " << lastB()
	<< " " << meanB()
	<< " " << stdDeviationB()
	<< " " << lastD()
	<< " " << meanD()
	<< " " << stdDeviationD()
	<< " " << lastObservable()
	<< " " << meanObservable()
	<< " " << stdDeviationObservable()
	<< " " << lastBetterObservable()
	<< " " << meanBetterObservable()
	<< " " << stdDeviationBetterObservable()
	<< " " << lastGeneratorOnBasis()
	<< " " << meanGeneratorOnBasis()
	<< " " << stdDeviationGeneratorOnBasis()   // 16
	//<< " " << -correlationB2()
	<< " " << (correlationB2() - meanB()) / (2*meanD())
	<< std::endl;*/

    out << time;
    for (size_t iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
      out << " " << meanB()(iOfFunction) * lastA(iOfFunction);
    for (size_t iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
      out << " " << lastA(iOfFunction);
    /*for (size_t iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
      out << " " << lastB()(iOfFunction);*/
    for (size_t iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
      out << " " << meanB()(iOfFunction);
    /*for (size_t iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
      for (size_t iOfFunction2 = 0; iOfFunction2 < nbOfFunctions(); iOfFunction2++)
	out << " " << lastD()(iOfFunction, iOfFunction2);*/
    for (size_t iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
      for (size_t iOfFunction2 = 0; iOfFunction2 < nbOfFunctions(); iOfFunction2++)
				out << " " << meanD()(iOfFunction, iOfFunction2);  //8-11

    out << " " << lastObservable()
      << " " << meanObservable()   //13
      << " " << stdDeviationObservable()
      << " " << lastBetterObservable()
      << " " << meanBetterObservable()
      << " " << stdDeviationBetterObservable();
    /*for (size_t iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
      out << " " << lastGeneratorOnBasis()(iOfFunction);*/
    for (size_t iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
      out<< " " << meanGeneratorOnBasis()(iOfFunction);

    for (size_t iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
      out << " " << meanB1()(iOfFunction);

    for (size_t iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
      out << " " << -correlationB2()(iOfFunction);  //22-23

    out << std::endl;

    /*for (size_t iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
    {
      out << " " <<
    }*/
  }

  void ControlVariate::postTreat(std::ofstream& out, double timeStep)
  {
      std::cout << "Post-treatment of the output" << std::endl;
    for (size_t iOfIteration = 0; iOfIteration < (size_t) historyObservable_.size(); iOfIteration++)
    {
      statsPostBetterObservable_.append(historyObservable_(iOfIteration) - lastA_.dot(historyGeneratorOnBasis_.row(iOfIteration)), iOfIteration);
      statsPostObservable_.append(historyObservable_(iOfIteration), iOfIteration);
			if (doOutput(iOfIteration))
			{
				out << iOfIteration * timeStep
						<< " " << statsPostObservable_.lastValue()
						<< " " << statsPostObservable_.mean()
						<< " " << statsPostObservable_.standardDeviation()
						<< " " << statsPostBetterObservable_.lastValue()
						<< " " << statsPostBetterObservable_.mean()
						<< " " << statsPostBetterObservable_.standardDeviation();
				for (size_t iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
					out << " " << lastA_(iOfFunction);
				out << std::endl;
			}
    }

    std::cout << "(D + Dt)/2 = " << std::endl << meanD() << std::endl;
    std::cout << "b = " << std::endl << meanB() << std::endl;
    std::cout << "a = " << std::endl << lastA_ << std::endl << std::endl;

    std::cout << "-<Linv j, j> = " << .5 * pow(statsPostObservable_.standardDeviation(), 2) << std::endl;
  }



  //############NoControlVariate###############"





  NoControlVariate::NoControlVariate(Input const& input, Potential& potential, size_t indexOfReplica):
    ControlVariate(input, potential, indexOfReplica, 1){}

  size_t NoControlVariate::nbOfFunctions() const
  {
    return 0;
  }

  size_t NoControlVariate::nbOfFunctionPairs() const
  {
    return 0;
  }

  bool NoControlVariate::isNone() const
  {
    return true;
  }

  double NoControlVariate::basisFunction(std::vector<Particle> const& /*configuration*/, size_t /*iOfFunction*/) const
  {
    return 0;
  }

    //double generatorOnBasisFunction(std::vector<Particle> const& configuration) const;
  double NoControlVariate::laplacienQ(std::vector<Particle> const& /*configuration*/, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    return 0;
  }

  Vector<double> NoControlVariate::gradientQ(std::vector<Particle> const& configuration, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    return Vector<double>(configuration[0].dimension());
  }

  double NoControlVariate::laplacienP(std::vector<Particle> const& /*configuration*/, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    return 0;
  }

  Vector<double> NoControlVariate::gradientP(std::vector<Particle> const& configuration, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    return Vector<double>(configuration[0].dimension());
  }

  void NoControlVariate::update(double observable, VectorXd& /*generatorOnBasisFunction*/, std::vector<Particle> const& /*configuration*/, size_t iOfIteration)
  {
    appendToObservable(observable, iOfIteration);
  }

  void NoControlVariate::postTreat(std::ofstream& out, double timeStep)
  {
      std::cout << "Post-treatment of the output" << std::endl;
    for (size_t iOfIteration = 0; iOfIteration < (size_t) historyObservable_.size(); iOfIteration++)
    {
      statsPostObservable_.append(historyObservable_(iOfIteration), iOfIteration);
      out << iOfIteration * timeStep
					<< " " << statsPostObservable_.lastValue()
					<< " " << statsPostObservable_.mean()
					<< " " << statsPostObservable_.standardDeviation()
					<< std::endl;
    }

    std::cout << "-<Linv j, j> = " << .5 * pow(statsPostObservable_.standardDeviation(), 2) << std::endl;
  }





  SinusControlVariate::SinusControlVariate(Input const& input, Potential& potential, size_t indexOfReplica):
    ControlVariate(input, potential, indexOfReplica, 1)
  {}

  double SinusControlVariate::basisFunction(std::vector<Particle> const& configuration, size_t /*iOfFunction*/) const
  {
    double q = configuration[0].position(0);
    return std::sin(2 * M_PI * q);
  }

  Vector<double> SinusControlVariate::gradientQ(std::vector<Particle> const& configuration, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    double q = configuration[0].position(0);
    Vector<double> grad(configuration[0].dimension());
    grad(0) = 2 * M_PI * std::cos(2 * M_PI * q);
    return grad;
  }

  double SinusControlVariate::laplacienQ(std::vector<Particle> const& configuration, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    double q = configuration[0].position(0);
    return - pow (2 * M_PI, 2) * std::sin(2 * M_PI * q);
  }


  double SinusControlVariate::laplacienP(std::vector<Particle> const& /*configuration*/, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    return 0;
  }

  Vector<double> SinusControlVariate::gradientP(std::vector<Particle> const& configuration, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    return Vector<double>(configuration[0].dimension());
  }




  CosControlVariate::CosControlVariate(Input const& input, Potential& potential, size_t indexOfReplica):
    ControlVariate(input, potential, indexOfReplica, 1)
  {}

  double CosControlVariate::basisFunction(std::vector<Particle> const& configuration, size_t /*iOfFunction*/) const
  {
    double q = configuration[0].position(0);
    return std::cos(2 * M_PI *q);
  }

  Vector<double> CosControlVariate::gradientQ(std::vector<Particle> const& configuration, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    double q = configuration[0].position(0);
    Vector<double> grad(configuration[0].dimension());
    grad(0) = - 2 * M_PI * std::sin(2 * M_PI * q);
    return grad;
  }

    double CosControlVariate::laplacienQ(std::vector<Particle> const& configuration, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    double q = configuration[0].position(0);
    return - pow (2 * M_PI, 2) * std::cos(2 * M_PI * q);
  }

  Vector<double> CosControlVariate::gradientP(std::vector<Particle> const& configuration, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    return Vector<double>(configuration[0].dimension());
  }

    double CosControlVariate::laplacienP(std::vector<Particle> const& /*configuration*/, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    return 0;
  }








  SinExpControlVariate::SinExpControlVariate(Input const& input, Potential& potential, size_t indexOfReplica):
    ControlVariate(input, potential, indexOfReplica, 1)
  {}

  double SinExpControlVariate::basisFunction(std::vector<Particle> const& configuration, size_t /*iOfFunction*/) const
  {
    Vector<double> q = configuration[0].position();
    double q0 = q(0);
    return std::sin(2 * M_PI * q0) * std::exp(potential(q)/2);
  }

  Vector<double> SinExpControlVariate::gradientQ(std::vector<Particle> const& configuration, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    Vector<double> q = configuration[0].position();
    double q0 = q(0);
    Vector<double> grad(configuration[0].dimension());
    grad(0) = (2 * M_PI * std::cos(2 * M_PI * q0)
	+ potentialDerivative(q)(0)/2 * std::sin(2 * M_PI * q0))
	* std::exp(potential(q)/2);
    return grad;
  }

    double SinExpControlVariate::laplacienQ(std::vector<Particle> const& configuration, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    Vector<double> q = configuration[0].position();
    double q0 = q(0);

    return (-pow(2 * M_PI, 2) * std::sin(2* M_PI * q0)
      +2 * M_PI * potentialDerivative(q)(0) *  std::cos(2 * M_PI * q0)
      + potentialLaplacian(q)/2 * std::sin(2 * M_PI * q0)
      + pow(potentialDerivative(q)(0) / 2, 2) * std::sin(2 * M_PI * q0))
      * std::exp(potential(q)/2);
  }

  Vector<double> SinExpControlVariate::gradientP(std::vector<Particle> const& configuration, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    return Vector<double>(configuration[0].dimension());
  }

  double SinExpControlVariate::laplacienP(std::vector<Particle> const& /*configuration*/, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    return 0;
  }




  CosExpControlVariate::CosExpControlVariate(Input const& input, Potential& potential, size_t indexOfReplica):
    ControlVariate(input, potential, indexOfReplica, 1)
  {}

  double CosExpControlVariate::basisFunction(std::vector<Particle> const& configuration, size_t /*iOfFunction*/) const
  {
    Vector<double> q = configuration[0].position();
    double q0 = q(0);
    return std::cos(2 * M_PI * q0) * std::exp(potential(q)/2);
  }

  Vector<double> CosExpControlVariate::gradientQ(std::vector<Particle> const& configuration, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    Vector<double> q = configuration[0].position();
    double q0 = q(0);
    Vector<double> grad(configuration[0].dimension());
    grad(0) = (- 2 * M_PI * std::sin(2 * M_PI * q0)
	+ potentialDerivative(q)(0)/2 * std::cos(2 * M_PI * q0))
	* std::exp(potential(q)/2);
    return grad;
  }

    double CosExpControlVariate::laplacienQ(std::vector<Particle> const& configuration, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    Vector<double> q = configuration[0].position();
    double q0 = q(0);

    return (-pow(2 * M_PI, 2) * std::cos(2* M_PI * q0)
      -2 * M_PI * potentialDerivative(q)(0) *  std::sin(2 * M_PI * q0)
      + potentialLaplacian(q)/2 * std::cos(2 * M_PI * q0)
      + pow(potentialDerivative(q)(0) / 2, 2) * std::cos(2 * M_PI * q0))
      * std::exp(potential(q)/2);
  }

  Vector<double> CosExpControlVariate::gradientP(std::vector<Particle> const& configuration, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    return Vector<double>(configuration[0].dimension());
  }

  double CosExpControlVariate::laplacienP(std::vector<Particle> const& /*configuration*/, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    return 0;
  }




  LangevinControlVariate::LangevinControlVariate(Input const& input, Potential& potential, size_t indexOfReplica):
    ControlVariate(input, potential, indexOfReplica, 1)
    {}

  double LangevinControlVariate::basisFunction(std::vector<Particle> const& configuration, size_t /*iOfFunction*/) const
  {
    //return pow(configuration[0].momentum(0), 2) * configuration[0].position(0);
    //return configuration[0].momentum(0);
    return configuration[0].kineticEnergy();
  }

  Vector<double> LangevinControlVariate::gradientQ(std::vector<Particle> const& configuration, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
     Vector<double> grad(configuration[0].dimension());
     //grad(0) = pow(configuration[0].momentum(0), 2);
     grad(0) = 0;
     return grad;
  }

    double LangevinControlVariate::laplacienQ(std::vector<Particle> const& /*configuration*/, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    assert(false);
    return 0;
  }

  Vector<double> LangevinControlVariate::gradientP(std::vector<Particle> const& configuration, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
     Vector<double> grad(configuration[0].dimension());
     grad(0) = configuration[0].momentum(0);
     //grad(0) = 1;
     return grad;
  }

  double LangevinControlVariate::laplacienP(std::vector<Particle> const& /*configuration*/, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    return 1;
    //return 2 * configuration[0].position(0);
  }



  SumEnergyControlVariate::SumEnergyControlVariate(Input const& input, Potential& potential, size_t indexOfReplica):
    ControlVariate(input, potential, indexOfReplica, 1), i0_(0)
  {}

  double SumEnergyControlVariate::basisFunction(std::vector<Particle> const& configuration, size_t /*iOfFunction*/) const
  {
    double result = 0;
    for (size_t iOfParticle = 1; iOfParticle < configuration.size(); iOfParticle++)
      result += configuration[iOfParticle].kineticEnergy() * (iOfParticle - i0_)
	    + (iOfParticle - i0_ - .5) * configuration[iOfParticle].potentialEnergy();
    //std::cout << result << std::endl;
    return result;
  }

  Vector<double> SumEnergyControlVariate::gradientQ(std::vector<Particle> const& configuration, size_t iOfParticle, size_t /*iOfFunction*/) const
  {
    return (iOfParticle - i0_ - .5) * configuration[iOfParticle].energyGrad();
  }

    double SumEnergyControlVariate::laplacienQ(std::vector<Particle> const& /*configuration*/, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    assert(false);
    return 0;
  }

  Vector<double> SumEnergyControlVariate::gradientP(std::vector<Particle> const& configuration, size_t iOfParticle, size_t /*iOfFunction*/) const
  {
    return (iOfParticle - i0_) * configuration[iOfParticle].momentum();
  }

  double SumEnergyControlVariate::laplacienP(std::vector<Particle> const& /*configuration*/, size_t iOfParticle, size_t /*iOfFunction*/) const
  {
    return iOfParticle - i0_;
  }




  EnergyControlVariate::EnergyControlVariate(Input const& input, Potential& potential, size_t indexOfReplica):
    ControlVariate(input, potential, indexOfReplica, 1), i0_(0)
  {}

  double EnergyControlVariate::basisFunction(std::vector<Particle> const& configuration, size_t /*iOfFunction*/) const
  {
    double result = 0;
    for (size_t iOfParticle = 0; iOfParticle < configuration.size(); iOfParticle++)
      result += configuration[iOfParticle].kineticEnergy()
	    + configuration[iOfParticle].potentialEnergy();
    //std::cout << result << std::endl;
    return result;
  }

  Vector<double> EnergyControlVariate::gradientQ(std::vector<Particle> const& configuration, size_t iOfParticle, size_t /*iOfFunction*/) const
  {
    return configuration[iOfParticle].energyGrad();
  }

    double EnergyControlVariate::laplacienQ(std::vector<Particle> const& /*configuration*/, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    assert(false);
    return 0;
  }

  Vector<double> EnergyControlVariate::gradientP(std::vector<Particle> const& configuration, size_t iOfParticle, size_t /*iOfFunction*/) const
  {
    return configuration[iOfParticle].momentum();
  }

  double EnergyControlVariate::laplacienP(std::vector<Particle> const& /*configuration*/, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    return 1;
  }



  LocalControlVariate::LocalControlVariate(Input const& input, Potential& potential, size_t indexOfReplica):
    ControlVariate(input, potential, indexOfReplica, 1)
  {}

  double LocalControlVariate::basisFunction(std::vector<Particle> const& configuration, size_t /*iOfFunction*/) const
  {
    //return configuration[0].kineticEnergy() + configuration[0].potentialEnergy();
    return configuration[0].energy();
  }

  Vector<double> LocalControlVariate::gradientQ(std::vector<Particle> const& configuration, size_t iOfParticle, size_t /*iOfFunction*/) const
  {
    Vector<double> result = Vector<double>(configuration[0].dimension());

    if (iOfParticle == 0)
      result(0) = std::sin(configuration[0].position(0));

    return result;
  }

    double LocalControlVariate::laplacienQ(std::vector<Particle> const& /*configuration*/, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    assert(false);
    return 0;
  }

  Vector<double> LocalControlVariate::gradientP(std::vector<Particle> const& configuration, size_t iOfParticle, size_t /*iOfFunction*/) const
  {
    Vector<double> result = Vector<double>(configuration[0].dimension());

    if (iOfParticle == 0)
      result(0) = configuration[0].momentum(0);

    return result;
  }

  double LocalControlVariate::laplacienP(std::vector<Particle> const& /*configuration*/, size_t iOfParticle, size_t /*iOfFunction*/) const
  {
    if (iOfParticle == 0)
      return 1;
    else
      return 0;
  }



  KineticControlVariate::KineticControlVariate(Input const& input, Potential& potential, size_t indexOfReplica):
    ControlVariate(input, potential, indexOfReplica, 1)
  {}

  double KineticControlVariate::basisFunction(std::vector<Particle> const& configuration, size_t /*iOfFunction*/) const
  {
    return 2 * configuration[0].momentum(0) * (std::sin(configuration[1].position(0) - configuration[0].position(0)) - std::sin(configuration[0].position(0)))
	  + 2 * configuration[0].kineticEnergy();
  }

  Vector<double> KineticControlVariate::gradientQ(std::vector<Particle> const& configuration, size_t iOfParticle, size_t /*iOfFunction*/) const
  {
    Vector<double> result(configuration[0].dimension());

    if (iOfParticle == 0)
      result(0) = -2 * configuration[0].momentum(0) * (std::cos(configuration[1].position(0) - configuration[0].position(0)) + std::cos(configuration[0].position(0)));
    else if (iOfParticle == 1)
      result(0) = 2 * configuration[0].momentum(0) * std::cos(configuration[1].position(0) - configuration[0].position(0));

    return result;
  }

    double KineticControlVariate::laplacienQ(std::vector<Particle> const& /*configuration*/, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    assert(false);
    return 0;
  }

  Vector<double> KineticControlVariate::gradientP(std::vector<Particle> const& configuration, size_t iOfParticle, size_t /*iOfFunction*/) const
  {
    Vector<double> result(configuration[0].dimension());

    if (iOfParticle == 0)
      result(0) = 2 * (std::sin(configuration[1].position(0) - configuration[0].position(0)) - std::sin(configuration[0].position(0)))
	  + 2 * configuration[0].momentum(0);

    return result;
  }

  double KineticControlVariate::laplacienP(std::vector<Particle> const& /*configuration*/, size_t iOfParticle, size_t /*iOfFunction*/) const
  {
    if (iOfParticle == 0)
      return 2;
    else
      return 0;
  }



  //#####TwoControlVariate#####

    TwoControlVariate::TwoControlVariate(Input const& input, Potential& potential, size_t indexOfReplica):
    ControlVariate(input, potential, indexOfReplica, 2)
  {}

  double TwoControlVariate::basisFunction(std::vector<Particle> const& configuration, size_t iOfFunction) const
  {
    /*if (iOfFunction == 0)
      return configuration[0].energy();
    else if (iOfFunction == 1)
      return configuration[0].momentum(0);
    else assert(false);*/

    if (iOfFunction == 0)
      return (configuration[0].momentum(0) - configuration[configuration.size()-1].momentum(0))
	  * (configuration[0].position(0) + configuration[configuration.size()-1].position(0));
    else if (iOfFunction == 1)
      return (configuration[0].momentum(0) + configuration[configuration.size()-1].momentum(0))
	  * (configuration[0].position(0) - configuration[configuration.size()-1].position(0));
    else assert(false);
  }

  Vector<double> TwoControlVariate::gradientQ(std::vector<Particle> const& configuration, size_t iOfParticle, size_t iOfFunction) const
  {
    Vector<double> result = Vector<double>(configuration[0].dimension());

    if (iOfFunction == 0)
      result(0) = - configuration[0].momentum(0) - configuration[configuration.size()-1].momentum(0);
    else if (iOfFunction == 1)
      {
      if (iOfParticle == 0)
	result(0) = configuration[0].momentum(0) + configuration[configuration.size()-1].momentum(0);
      else if (iOfParticle == 1)
	result(0) = - configuration[0].momentum(0) - configuration[configuration.size()-1].momentum(0);
    }
      return result;
  }

    double TwoControlVariate::laplacienQ(std::vector<Particle> const& /*configuration*/, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    assert(false);
    return 0;
  }

  Vector<double> TwoControlVariate::gradientP(std::vector<Particle> const& configuration, size_t iOfParticle, size_t iOfFunction) const
  {
    Vector<double> result(configuration[0].dimension());

    if (iOfFunction == 0)
    {
      if (iOfParticle == 0)
	result(0) = configuration[0].position(0) + configuration[configuration.size()-1].position(0);
      else if (iOfParticle == 1)
	result(0) = - configuration[0].position(0) - configuration[configuration.size()-1].position(0);
    }
    else if (iOfFunction == 1)
      result(0) = - configuration[0].position(0) - configuration[configuration.size()-1].position(0);


    return result;
  }

  double TwoControlVariate::laplacienP(std::vector<Particle> const& /*configuration*/, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    return 0;
  }




   //#####IOControlVariate#####

  IOControlVariate::IOControlVariate(Input const& input, Potential& potential, size_t indexOfReplica):
    ControlVariate(input, potential, indexOfReplica, 2)
  {}

  double IOControlVariate::basisFunction(std::vector<Particle> const& configuration, size_t iOfFunction) const
  {
    /*if (iOfFunction == 0)
      return configuration[0].momentum(0) * configuration[1].momentum(0);
    else if (iOfFunction == 1)
      return configuration[0].momentum(0) * configuration[2].momentum(0);
    else assert(false);*/

    if (iOfFunction == 0)
      return configuration[0].energy();
    else if (iOfFunction == 1)
      return -configuration[0].momentum(0) * configuration[0].energyGrad(0);
    else assert(false);
  }

  Vector<double> IOControlVariate::gradientQ(std::vector<Particle> const& configuration, size_t iOfParticle, size_t iOfFunction) const
  {
    Vector<double> result = Vector<double>(configuration[0].dimension());

    if (iOfFunction == 0 && iOfParticle == 0)
      result(0) = std::sin(configuration[0].position(0));
    else if (iOfFunction == 1)
    {
      if (iOfParticle == 0)
	result(0) = configuration[0].momentum(0) * std::cos(configuration[1].position(0) - configuration[0].position(0));
      else if (iOfParticle == 1)
	result(0) = -configuration[0].momentum(0) * std::cos(configuration[1].position(0) - configuration[0].position(0));
    }
    return result;
  }

    double IOControlVariate::laplacienQ(std::vector<Particle> const& /*configuration*/, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    assert(false);
    return 0;
  }

  Vector<double> IOControlVariate::gradientP(std::vector<Particle> const& configuration, size_t iOfParticle, size_t iOfFunction) const
  {
    Vector<double> result(configuration[0].dimension());
    /*if (iOfFunction == 0)
    {
      if (iOfParticle == 0)
	return configuration[1].momentum();
      else if (iOfParticle == 1)
	return configuration[0].momentum();
    }
    else if (iOfFunction == 1)
    {
      if (iOfParticle == 0)
	return configuration[2].momentum();
      else if (iOfParticle == 2)
	return configuration[0].momentum();
    }
    else assert(false);*/

    if (iOfFunction == 0 && iOfParticle == 0)
      result(0) = configuration[0].momentum(0);
    else if (iOfFunction == 1 && iOfParticle == 0)
      result(0) = - configuration[0].energyGrad(0);

    return result;
  }

  double IOControlVariate::laplacienP(std::vector<Particle> const& /*configuration*/, size_t iOfParticle, size_t iOfFunction) const
  {
    if (iOfFunction == 0 && iOfParticle == 0)
      return 1;
    else
      return 0;
  }

}
