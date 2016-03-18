#include "controlVariate.hpp"

using std::cout; 
using std::endl; 
using std::sin;
using std::cos;
using std::exp;
using std::ofstream;

namespace simol 
{
  
  ControlVariate* createControlVariate(Input const& input, Potential& potential, Galerkin* galerkin)
  {
    if (input.controlVariateName() == "None")
      return new NoControlVariate(input, potential);
    if (input.controlVariateName() == "Sinus")
      return new SinusControlVariate(input, potential);
    if (input.controlVariateName() == "Cosine")
      return new CosControlVariate(input, potential);
    else if (input.controlVariateName() == "SinExp")
      return new SinExpControlVariate(input, potential);
    else if (input.controlVariateName() == "CosExp")
      return new CosExpControlVariate(input, potential);
    else if (input.controlVariateName() == "Langevin")
      return new LangevinControlVariate(input, potential);
    else if (input.controlVariateName() == "SumEnergy")
      return new SumEnergyControlVariate(input, potential);
    else if (input.controlVariateName() == "Energy")
      return new EnergyControlVariate(input, potential);
    else if (input.controlVariateName() == "Local")
      return new LocalControlVariate(input, potential);
    else if (input.controlVariateName() == "Kinetic")
      return new KineticControlVariate(input, potential);
    else if (input.controlVariateName() == "Two")
      return new TwoControlVariate(input, potential);
		else if (input.controlVariateName() == "ExpFourierHermite")
			return new ExpFourierHermiteControlVariate(input, potential, galerkin);
    else
      cout << input.controlVariateName() << " is not a valid control variate !" << std::endl;
    return 0;
  }

  
  ControlVariate::ControlVariate(Input const& input, Potential& potential, size_t nbOfFunctions):
    decorrelationNbOfIterations_(input.decorrelationNbOfIterations()),
    decorrelationTime_(input.decorrelationTime()),
    nbOfFunctions_(nbOfFunctions),
    nbOfFunctionPairs_(pow(nbOfFunctions_, 2)),
    periodNbOfIterations_(input.outputPeriodNbOfIterations()),
    nbOfAutocoPts_(input.nbOfAutocoPts()),
    statsObservable_(decorrelationNbOfIterations(), decorrelationTime(), nbOfAutocoPts()),
    statsBetterObservable_(decorrelationNbOfIterations(), decorrelationTime(), nbOfAutocoPts()),
    statsGeneratorOnBasis_(nbOfFunctions_),
    statsB1_(nbOfFunctions_),
    statsB2_(decorrelationNbOfIterations(), decorrelationTime(), nbOfAutocoPts(), nbOfFunctions_),
    statsD_(nbOfFunctions_, nbOfFunctions_),
    lastA_(nbOfFunctions_),
    statsPostObservable_(decorrelationNbOfIterations(), decorrelationTime(), nbOfAutocoPts()),
    statsPostBetterObservable_(decorrelationNbOfIterations(), decorrelationTime(), nbOfAutocoPts()),
    historyObservable_(input.nbOfIterations()),
    historyGeneratorOnBasis_(input.nbOfIterations(), nbOfFunctions_),
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
  
  bool ControlVariate::doOutput(size_t iOfIteration) const
  {
    return (periodNbOfIterations_ > 0 && iOfIteration % periodNbOfIterations_ == 0);
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

  
  
  size_t ControlVariate::decorrelationNbOfIterations() const
  {
    return decorrelationNbOfIterations_;
  }
  
  double ControlVariate::decorrelationTime() const
  {
    return decorrelationTime_;
  }
  
  int const& ControlVariate::nbOfAutocoPts() const
  {
		return nbOfAutocoPts_;
	}
	
	int ControlVariate::nbOfFourier() const
	{return 0;}
  
  int ControlVariate::nbOfHermite() const
	{return 0;}
	
	
	
	
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
    //Eigen::Matrix<AutocorrelationStats<double>, Eigen::Dynamic, Eigen::Dynamic> A(nbOfFunctions_, nbOfFunctions_, AutocorrelationStats<double>(decorrelationNbOfIterations(), decorrelationTime()));
    
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
    cout << B.transpose() * Dinv << endl;*/
    //cout << (statsB_.meanMat() - statsB2_.integratedAutocorrelationMat()).transpose().size() << "  " << statsD_.meanMat().inverse().size() << endl;
    //lastA_ = (statsB_.meanMat() - statsB2_.integratedAutocorrelationMat()).transpose() * statsD_.meanMat().llt().solve(generatorOnBasisFunction);
    //cout << lastA_.size() << "   " << generatorOnBasisFunction.size() << endl;
    
    if (statsD_.meanMat().determinant() != 0)
    {
      //lastA_ = - .5 * statsD_.meanMat().llt().solve(meanB());
			//lastA_ = VectorXd(1,1);
			lastA_.fill(1);
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

  
  void ControlVariate::update(double observable, VectorXd& generatorOnBasisFunction, vector<Particle> const& configuration, size_t iOfIteration)
  {
		//cout << "ControlVariate::update" << endl;
    VectorXd valueBasisFunction(nbOfFunctions_);
    for (size_t iOfFunction =0; iOfFunction < nbOfFunctions_; iOfFunction++)
      valueBasisFunction(iOfFunction) = basisFunction(configuration, iOfFunction);
    appendToB1(observable, valueBasisFunction);
    appendToB2(observable, generatorOnBasisFunction, iOfIteration);
    appendToD(generatorOnBasisFunction, valueBasisFunction);
    appendToObservable(observable, iOfIteration);
    appendToBetterObservable(observable, generatorOnBasisFunction, iOfIteration);
		//cout << "end ControlVariate::update" << endl;
  }
    
    void ControlVariate::display(std::ofstream& out, double time) const
  {    
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
    for (size_t iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
      out << " " << lastGeneratorOnBasis()(iOfFunction);
    for (size_t iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
      out<< " " << meanGeneratorOnBasis()(iOfFunction);
    
    for (size_t iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
      out << " " << meanB1()(iOfFunction);
    
    for (size_t iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
      out << " " << -correlationB2()(iOfFunction);  //22-23
    
    out << endl;
	
    /*for (size_t iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
    {
      out << " " << 
    }*/
  }
  
  void ControlVariate::postTreat(std::ofstream& out, double timeStep)
  {
    cout << "Post-treatment of the output" << endl;
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
				out << endl;
			}
    }
    
    cout << "(D + Dt)/2 = " << endl << meanD() << endl;
    cout << "b = " << endl << meanB() << endl;
    cout << "a = " << endl << lastA_ << endl << endl;
    
    cout << "-<Linv j, j> = " << .5 * pow(statsPostObservable_.standardDeviation(), 2) << endl;
  }
  
 
  
  //############NoControlVariate###############"
  
  
  
  
  
  NoControlVariate::NoControlVariate(Input const& input, Potential& potential):
    ControlVariate(input, potential, 1){}
    
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
    
  double NoControlVariate::basisFunction(vector<Particle> const& /*configuration*/, size_t /*iOfFunction*/) const
  {
    return 0;
  }
  
    //double generatorOnBasisFunction(vector<Particle> const& configuration) const;
  double NoControlVariate::laplacianQ(vector<Particle> const& /*configuration*/, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    return 0;
  }
  
  Vector<double> NoControlVariate::gradientQ(vector<Particle> const& configuration, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    return Vector<double>(configuration[0].dimension());
  }
  
  double NoControlVariate::laplacianP(vector<Particle> const& /*configuration*/, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    return 0;
  }
    
  Vector<double> NoControlVariate::gradientP(vector<Particle> const& configuration, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    return Vector<double>(configuration[0].dimension());
  }
  
  void NoControlVariate::update(double observable, VectorXd& /*generatorOnBasisFunction*/, vector<Particle> const& /*configuration*/, size_t iOfIteration)
  {
    appendToObservable(observable, iOfIteration);
  }

  void NoControlVariate::postTreat(std::ofstream& out, double timeStep)
  {
    cout << "Post-treatment of the output" << endl;
    for (size_t iOfIteration = 0; iOfIteration < (size_t) historyObservable_.size(); iOfIteration++)
    {
      statsPostObservable_.append(historyObservable_(iOfIteration), iOfIteration);
			if (doOutput(iOfIteration))
			{
				out << iOfIteration * timeStep
						<< " " << statsPostObservable_.lastValue() 
						<< " " << statsPostObservable_.mean() 
						<< " " << statsPostObservable_.standardDeviation() 
						<< endl;
			}
    }
    
    cout << "-<Linv j, j> = " << .5 * pow(statsPostObservable_.standardDeviation(), 2) << endl;
  } 
    
    
  
  
  
  SinusControlVariate::SinusControlVariate(Input const& input, Potential& potential):
    ControlVariate(input, potential, 1)
  {}
  
  double SinusControlVariate::basisFunction(vector<Particle> const& configuration, size_t /*iOfFunction*/) const
  {
    double q = configuration[0].position(0);
    return sin(2 * M_PI * q);
  }
   
  Vector<double> SinusControlVariate::gradientQ(vector<Particle> const& configuration, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    double q = configuration[0].position(0);
    Vector<double> grad(configuration[0].dimension());
    grad(0) = 2 * M_PI * cos(2 * M_PI * q);
    return grad;
  }
  
  double SinusControlVariate::laplacianQ(vector<Particle> const& configuration, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    double q = configuration[0].position(0);
    return - pow (2 * M_PI, 2) * sin(2 * M_PI * q);
  }

  
  double SinusControlVariate::laplacianP(vector<Particle> const& /*configuration*/, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    return 0;
  }
  
  Vector<double> SinusControlVariate::gradientP(vector<Particle> const& configuration, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    return Vector<double>(configuration[0].dimension());
  }
  
  
  
  
  CosControlVariate::CosControlVariate(Input const& input, Potential& potential):
    ControlVariate(input, potential, 1)
  {}
  
  double CosControlVariate::basisFunction(vector<Particle> const& configuration, size_t /*iOfFunction*/) const
  {
    double q = configuration[0].position(0);
    return cos(2 * M_PI *q);
  }
  
  Vector<double> CosControlVariate::gradientQ(vector<Particle> const& configuration, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    double q = configuration[0].position(0);
    Vector<double> grad(configuration[0].dimension());
    grad(0) = - 2 * M_PI * sin(2 * M_PI * q);
    return grad;
  }
  
    double CosControlVariate::laplacianQ(vector<Particle> const& configuration, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    double q = configuration[0].position(0);
    return - pow (2 * M_PI, 2) * cos(2 * M_PI * q);
  }
  
  Vector<double> CosControlVariate::gradientP(vector<Particle> const& configuration, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    return Vector<double>(configuration[0].dimension());
  }
  
    double CosControlVariate::laplacianP(vector<Particle> const& /*configuration*/, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    return 0;
  }
  
  
  
  
  
  
  
  
  SinExpControlVariate::SinExpControlVariate(Input const& input, Potential& potential):
    ControlVariate(input, potential, 1)
  {}
  
  double SinExpControlVariate::basisFunction(vector<Particle> const& configuration, size_t /*iOfFunction*/) const
  {
    Vector<double> q = configuration[0].position();
    double q0 = q(0);
    return sin(2 * M_PI * q0) * exp(potential(q)/2);
  }
  
  Vector<double> SinExpControlVariate::gradientQ(vector<Particle> const& configuration, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    Vector<double> q = configuration[0].position();
    double q0 = q(0);
    Vector<double> grad(configuration[0].dimension());
    grad(0) = (2 * M_PI * cos(2 * M_PI * q0)
	+ potentialDerivative(q)(0)/2 * sin(2 * M_PI * q0))
	* exp(potential(q)/2);
    return grad;
  }
  
    double SinExpControlVariate::laplacianQ(vector<Particle> const& configuration, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    Vector<double> q = configuration[0].position();
    double q0 = q(0);
    
    return (-pow(2 * M_PI, 2) * sin(2* M_PI * q0)
      +2 * M_PI * potentialDerivative(q)(0) *  cos(2 * M_PI * q0)
      + potentialLaplacian(q)/2 * sin(2 * M_PI * q0)
      + pow(potentialDerivative(q)(0) / 2, 2) * sin(2 * M_PI * q0))
      * exp(potential(q)/2);
  }
  
  Vector<double> SinExpControlVariate::gradientP(vector<Particle> const& configuration, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    return Vector<double>(configuration[0].dimension());
  }
    
  double SinExpControlVariate::laplacianP(vector<Particle> const& /*configuration*/, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    return 0;
  }
  
  
  
  
  CosExpControlVariate::CosExpControlVariate(Input const& input, Potential& potential):
    ControlVariate(input, potential, 1)
  {}
  
  double CosExpControlVariate::basisFunction(vector<Particle> const& configuration, size_t /*iOfFunction*/) const
  {
    Vector<double> q = configuration[0].position();
    double q0 = q(0);
    return cos(2 * M_PI * q0) * exp(potential(q)/2);
  }
  
  Vector<double> CosExpControlVariate::gradientQ(vector<Particle> const& configuration, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    Vector<double> q = configuration[0].position();
    double q0 = q(0);
    Vector<double> grad(configuration[0].dimension());
    grad(0) = (- 2 * M_PI * sin(2 * M_PI * q0)
	+ potentialDerivative(q)(0)/2 * cos(2 * M_PI * q0))
	* exp(potential(q)/2);
    return grad;
  }
  
    double CosExpControlVariate::laplacianQ(vector<Particle> const& configuration, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    Vector<double> q = configuration[0].position();
    double q0 = q(0);
    
    return (-pow(2 * M_PI, 2) * cos(2* M_PI * q0)
      -2 * M_PI * potentialDerivative(q)(0) *  sin(2 * M_PI * q0)
      + potentialLaplacian(q)/2 * cos(2 * M_PI * q0)
      + pow(potentialDerivative(q)(0) / 2, 2) * cos(2 * M_PI * q0))
      * exp(potential(q)/2);
  }
  
  Vector<double> CosExpControlVariate::gradientP(vector<Particle> const& configuration, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    return Vector<double>(configuration[0].dimension());
  }
    
  double CosExpControlVariate::laplacianP(vector<Particle> const& /*configuration*/, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    return 0;
  }
  
  
  
  
  LangevinControlVariate::LangevinControlVariate(Input const& input, Potential& potential):
    ControlVariate(input, potential, 1)
    {}
  
  double LangevinControlVariate::basisFunction(vector<Particle> const& configuration, size_t /*iOfFunction*/) const
  {
    //return pow(configuration[0].momentum(0), 2) * configuration[0].position(0);
    //return configuration[0].momentum(0);
    return configuration[0].kineticEnergy();
  }
  
  Vector<double> LangevinControlVariate::gradientQ(vector<Particle> const& configuration, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
     Vector<double> grad(configuration[0].dimension());
     //grad(0) = pow(configuration[0].momentum(0), 2);
     grad(0) = 0;
     return grad;
  }
  
    double LangevinControlVariate::laplacianQ(vector<Particle> const& /*configuration*/, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    assert(false);
    return 0;
  }
  
  Vector<double> LangevinControlVariate::gradientP(vector<Particle> const& configuration, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
     Vector<double> grad(configuration[0].dimension());
     grad(0) = configuration[0].momentum(0);
     //grad(0) = 1;
     return grad;
  }
    
  double LangevinControlVariate::laplacianP(vector<Particle> const& /*configuration*/, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    return 1;
    //return 2 * configuration[0].position(0);
  }
  
  
  
  SumEnergyControlVariate::SumEnergyControlVariate(Input const& input, Potential& potential):
    ControlVariate(input, potential, 1), i0_(0)
  {}
  
  double SumEnergyControlVariate::basisFunction(vector<Particle> const& configuration, size_t /*iOfFunction*/) const
  {
    double result = 0;
    for (size_t iOfParticle = 1; iOfParticle < configuration.size(); iOfParticle++)
      result += configuration[iOfParticle].kineticEnergy() * (iOfParticle - i0_)
	    + (iOfParticle - i0_ - .5) * configuration[iOfParticle].potentialEnergy();
    //cout << result << endl;
    return result;
  }
  
  Vector<double> SumEnergyControlVariate::gradientQ(vector<Particle> const& configuration, size_t iOfParticle, size_t /*iOfFunction*/) const
  {
    return (iOfParticle - i0_ - .5) * configuration[iOfParticle].energyGrad();
  }
  
    double SumEnergyControlVariate::laplacianQ(vector<Particle> const& /*configuration*/, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    assert(false);
    return 0;
  }
  
  Vector<double> SumEnergyControlVariate::gradientP(vector<Particle> const& configuration, size_t iOfParticle, size_t /*iOfFunction*/) const
  {
    return (iOfParticle - i0_) * configuration[iOfParticle].momentum();
  }
    
  double SumEnergyControlVariate::laplacianP(vector<Particle> const& /*configuration*/, size_t iOfParticle, size_t /*iOfFunction*/) const
  {
    return iOfParticle - i0_;
  }
  
  
  
  
  EnergyControlVariate::EnergyControlVariate(Input const& input, Potential& potential):
    ControlVariate(input, potential, 1), i0_(0)
  {}
  
  double EnergyControlVariate::basisFunction(vector<Particle> const& configuration, size_t /*iOfFunction*/) const
  {
    double result = 0;
    for (size_t iOfParticle = 0; iOfParticle < configuration.size(); iOfParticle++)
      result += configuration[iOfParticle].kineticEnergy()
	    + configuration[iOfParticle].potentialEnergy();
    //cout << result << endl;
    return result;
  }
  
  Vector<double> EnergyControlVariate::gradientQ(vector<Particle> const& configuration, size_t iOfParticle, size_t /*iOfFunction*/) const
  {
    return configuration[iOfParticle].energyGrad();
  }
  
    double EnergyControlVariate::laplacianQ(vector<Particle> const& /*configuration*/, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    assert(false);
    return 0;
  }
  
  Vector<double> EnergyControlVariate::gradientP(vector<Particle> const& configuration, size_t iOfParticle, size_t /*iOfFunction*/) const
  {
    return configuration[iOfParticle].momentum();
  }
    
  double EnergyControlVariate::laplacianP(vector<Particle> const& /*configuration*/, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    return 1;
  }
  
  
  
  LocalControlVariate::LocalControlVariate(Input const& input, Potential& potential):
    ControlVariate(input, potential, 1)
  {}
  
  double LocalControlVariate::basisFunction(vector<Particle> const& configuration, size_t /*iOfFunction*/) const
  {
    //return configuration[0].kineticEnergy() + configuration[0].potentialEnergy();
    return configuration[0].energy(); 
  }
  
  Vector<double> LocalControlVariate::gradientQ(vector<Particle> const& configuration, size_t iOfParticle, size_t /*iOfFunction*/) const
  {
    Vector<double> result = Vector<double>(configuration[0].dimension());
    
    if (iOfParticle == 0)
      result(0) = sin(configuration[0].position(0));
    
    return result;
  }
  
    double LocalControlVariate::laplacianQ(vector<Particle> const& /*configuration*/, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    assert(false);
    return 0;
  }
  
  Vector<double> LocalControlVariate::gradientP(vector<Particle> const& configuration, size_t iOfParticle, size_t /*iOfFunction*/) const
  {
    Vector<double> result = Vector<double>(configuration[0].dimension());
    
    if (iOfParticle == 0)
      result(0) = configuration[0].momentum(0);
    
    return result;
  }
    
  double LocalControlVariate::laplacianP(vector<Particle> const& /*configuration*/, size_t iOfParticle, size_t /*iOfFunction*/) const
  {
    if (iOfParticle == 0)
      return 1;
    else 
      return 0;
  }
  
  
  
  KineticControlVariate::KineticControlVariate(Input const& input, Potential& potential):
    ControlVariate(input, potential, 1)
  {}
  
  double KineticControlVariate::basisFunction(vector<Particle> const& configuration, size_t /*iOfFunction*/) const
  {
    return 2 * configuration[0].momentum(0) * (sin(configuration[1].position(0) - configuration[0].position(0)) - sin(configuration[0].position(0)))
	  + 2 * configuration[0].kineticEnergy();
  }
  
  Vector<double> KineticControlVariate::gradientQ(vector<Particle> const& configuration, size_t iOfParticle, size_t /*iOfFunction*/) const
  {
    Vector<double> result(configuration[0].dimension());
    
    if (iOfParticle == 0)
      result(0) = -2 * configuration[0].momentum(0) * (cos(configuration[1].position(0) - configuration[0].position(0)) + cos(configuration[0].position(0)));
    else if (iOfParticle == 1)
      result(0) = 2 * configuration[0].momentum(0) * cos(configuration[1].position(0) - configuration[0].position(0));
    
    return result;
  }
  
    double KineticControlVariate::laplacianQ(vector<Particle> const& /*configuration*/, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    assert(false);
    return 0;
  }
  
  Vector<double> KineticControlVariate::gradientP(vector<Particle> const& configuration, size_t iOfParticle, size_t /*iOfFunction*/) const
  {
    Vector<double> result(configuration[0].dimension());
    
    if (iOfParticle == 0)
      result(0) = 2 * (sin(configuration[1].position(0) - configuration[0].position(0)) - sin(configuration[0].position(0)))
	  + 2 * configuration[0].momentum(0);
	  
    return result;
  }
    
  double KineticControlVariate::laplacianP(vector<Particle> const& /*configuration*/, size_t iOfParticle, size_t /*iOfFunction*/) const
  {
    if (iOfParticle == 0)
      return 2;
    else
      return 0;
  }
  
  
  
  //#####TwoControlVariate#####
  
    TwoControlVariate::TwoControlVariate(Input const& input, Potential& potential):
    ControlVariate(input, potential, 2)
  {}
  
  double TwoControlVariate::basisFunction(vector<Particle> const& configuration, size_t iOfFunction) const
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
		else
			throw std::invalid_argument("iOfFunction must 0 or 1");
  }
  
  Vector<double> TwoControlVariate::gradientQ(vector<Particle> const& configuration, size_t iOfParticle, size_t iOfFunction) const
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
  
    double TwoControlVariate::laplacianQ(vector<Particle> const& /*configuration*/, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    assert(false);
    return 0;
  }
  
  Vector<double> TwoControlVariate::gradientP(vector<Particle> const& configuration, size_t iOfParticle, size_t iOfFunction) const
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
    
  double TwoControlVariate::laplacianP(vector<Particle> const& /*configuration*/, size_t /*iOfParticle*/, size_t /*iOfFunction*/) const
  {
    return 0;
  }
  
  
  //#####BasisControlVariate#####
  
  BasisControlVariate::BasisControlVariate(Input const& input, Potential& potential, Galerkin* galerkin):
    ControlVariate(input, potential, 1),
    coeffsVec_(input.nbOfFourier() * input.nbOfHermite(), 1)
  {
		if (galerkin)
		{
			cout << "CVcoeffs from Galerkin !" << endl;
			coeffsVec_ = galerkin->CVcoeffs();
		}
		else
		{
			cout << "CVcoeffs from file !" << endl;
			DMat temp;
			temp.load(input.outputFolderName() + input.controlVariateCoeffsPath());
			assert(temp.n_rows == coeffsVec_.n_rows && temp.n_cols == 1);
			coeffsVec_ = arma::conv_to<SMat>::from(temp);
		}
		//cout << "coeffsVec_ : " << endl;
		//cout << coeffsVec_ << endl;
	}
  
  double BasisControlVariate::basisFunction(vector<Particle> const& configuration, size_t iOfFunction) const
  {
		assert(iOfFunction == 0);
		double result = 0;
		
		for (SMat::iterator it = coeffsVec_.begin_col(0); it != coeffsVec_.end_col(0); ++it)
		{
			int iOfCoeff = it.row();
			double valOfCoeff = *it;
			result += valOfCoeff * basis_->value(configuration, iOfCoeff);
		}
		
		return result;
  }
  
  Vector<double> BasisControlVariate::gradientQ(vector<Particle> const& configuration, size_t iOfParticle, size_t iOfFunction) const
  {
		//cout << "BasisControlVariate::gradientQ" << endl;
    Vector<double> result(1, 0);
    assert(iOfFunction == 0);
		//cout << "gradientQ" << endl;
		for (SMat::iterator it = coeffsVec_.begin_col(0); it != coeffsVec_.end_col(0); ++it)
		{
			int iOfCoeff = it.row();
			double valOfCoeff = *it;
			//cout << "+ " << valOfCoeff << " X ";
			result += valOfCoeff * basis_->gradientQ(configuration, iOfParticle, iOfCoeff);
		}
		//cout << endl;
    //cout << "end BasisControlVariate::gradientQ" << endl;
    //cout << "gradientQ = " << result << endl;
    return result;
  }
  
  double BasisControlVariate::laplacianQ(vector<Particle> const& configuration, size_t iOfParticle, size_t iOfFunction) const
  {
    double result = 0;
    assert(iOfFunction == 0);
		
		for (SMat::iterator it = coeffsVec_.begin_col(0); it != coeffsVec_.end_col(0); ++it)
		{
			int iOfCoeff = it.row();
			double valOfCoeff = *it;
			result += valOfCoeff * basis_->laplacianQ(configuration, iOfParticle, iOfCoeff);
		}
    //cout << "laplacianQ = " << result << endl;
    return result;
  }
  
  Vector<double> BasisControlVariate::gradientP(vector<Particle> const& configuration, size_t iOfParticle, size_t iOfFunction) const
  {
		Vector<double> result(1, 0);
    assert(iOfFunction == 0);
		//cout << "gradientP" << endl;
		for (SMat::iterator it = coeffsVec_.begin_col(0); it != coeffsVec_.end_col(0); ++it)
		{
			int iOfCoeff = it.row();
			double valOfCoeff = *it;
			//cout << "+ " << valOfCoeff << " X ";
			result += valOfCoeff * basis_->gradientP(configuration, iOfParticle, iOfCoeff);
		}
    //cout << endl;
		//cout << "gradientP = " << result << endl;
    return result;
  }
    
  double BasisControlVariate::laplacianP(vector<Particle> const& configuration, size_t iOfParticle, size_t iOfFunction) const
  {
    double result = 0;
    assert(iOfFunction == 0);
		
		for (SMat::iterator it = coeffsVec_.begin_col(0); it != coeffsVec_.end_col(0); ++it)
		{
			int iOfCoeff = it.row();
			double valOfCoeff = *it;
			result += valOfCoeff * basis_->laplacianP(configuration, iOfParticle, iOfCoeff);
		}
    //cout << "laplacianP = " << result << endl;
    return result;
  }
  
  
  
  ExpFourierHermiteControlVariate::ExpFourierHermiteControlVariate(Input const& input, Potential& potential, Galerkin* galerkin):
		BasisControlVariate(input, potential, galerkin)
	{		
		basis_ = new ExpFourierHermiteBasis(input, potential);
		nbQ_ = 100;  // rafinement du maillage de l'output map
		nbP_ = 100;  // rafinement du maillage de l'output map
		deltaQ_ = 2*M_PI/nbQ_;
		pMax_ = 4;
		deltaP_ = 2*pMax_/nbP_;
	}
	
	int ExpFourierHermiteControlVariate::nbOfFourier() const
	{return basis_->nbOfElts(0);}
  
  int ExpFourierHermiteControlVariate::nbOfHermite() const
	{return basis_->nbOfElts(1);}
	
	void ExpFourierHermiteControlVariate::displayMap(ofstream& out) const
	{
		cout << "displayMap" << endl;
		vector<Particle> conf(1);
		conf[0] = Particle(0,0,0);
		out << nbP_ + 1 << " ";
		for (double p = -pMax_; p<pMax_; p += deltaP_)
			out << p << " ";
		out << endl;
		for (double q = -M_PI; q < M_PI; q += deltaQ_)
		{
			out << q << " ";
			for (double p = -pMax_; p<pMax_; p += deltaP_)
			{
				conf[0].position() = Vector<double>(1, q);
				conf[0].momentum() = Vector<double>(1, p);
				out << basisFunction(conf, 0) << " ";
			}
			out << endl;
		}
		cout << "end displayMap" << endl;
	}
	
	void ExpFourierHermiteControlVariate::displayGradQMap(ofstream& out) const
		{
		cout << "displayGradQMap" << endl;
		vector<Particle> conf(1);
		conf[0] = Particle(0,0,0);
		out << nbP_ + 1 << " ";
		for (double p = -pMax_; p<pMax_; p += deltaP_)
			out << p << " ";
		out << endl;
		for (double q = -M_PI; q < M_PI; q += deltaQ_)
		{
			out << q << " ";
			for (double p = -pMax_; p<pMax_; p += deltaP_)
			{
				conf[0].position() = Vector<double>(1, q);
				conf[0].momentum() = Vector<double>(1, p);
				out << gradientQ(conf, 0) << " ";
			}
			out << endl;
		}
		cout << "end displayMap" << endl;
	}
	
		void ExpFourierHermiteControlVariate::displayGradPMap(ofstream& out) const
		{
		cout << "displayGradPMap" << endl;
		vector<Particle> conf(1);
		conf[0] = Particle(0,0,0);
		out << nbP_ + 1 << " ";
		for (double p = -pMax_; p<pMax_; p += deltaP_)
			out << p << " ";
		out << endl;
		for (double q = -M_PI; q < M_PI; q += deltaQ_)
		{
			out << q << " ";
			for (double p = -pMax_; p<pMax_; p += deltaP_)
			{
				conf[0].position() = Vector<double>(1, q);
				conf[0].momentum() = Vector<double>(1, p);
				out << gradientP(conf, 0) << " ";
			}
			out << endl;
		}
		cout << "end displayMap" << endl;
	}
	
}










