#include "simol/statphys/controlVariate/ControlVariate.hpp"

using std::cout;
using std::endl;
using std::sin;
using std::cos;
using std::exp;
using std::ofstream;

namespace simol
{

  Observable* createControlVariate(Input const& input, const string& outPath, Potential& potential, Galerkin* galerkin)
  {
    if (input.controlVariateName() == "None")
      return new Observable(input, outPath);
    if (input.controlVariateName() == "Sinus")
      return new SinusControlVariate(input, outPath, potential);
    if (input.controlVariateName() == "Cosine")
      return new CosControlVariate(input, outPath, potential);
    else if (input.controlVariateName() == "SinExp")
      return new SinExpControlVariate(input, outPath, potential);
    else if (input.controlVariateName() == "CosExp")
      return new CosExpControlVariate(input, outPath, potential);
    else if (input.controlVariateName() == "Langevin")
      return new LangevinControlVariate(input, outPath, potential);
    else if (input.controlVariateName() == "SumEnergy")
      return new SumEnergyControlVariate(input, outPath, potential);
    else if (input.controlVariateName() == "Energy")
      return new EnergyControlVariate(input, outPath, potential);
    else if (input.controlVariateName() == "Local")
      return new LocalControlVariate(input, outPath, potential);
    else if (input.controlVariateName() == "Kinetic")
      return new KineticControlVariate(input, outPath, potential);
    else if (input.controlVariateName() == "Two")
      return new TwoControlVariate(input, outPath, potential);
    else if (input.controlVariateName() == "ExpFourierHermite")
      return new ExpFourierHermiteControlVariate(input, outPath, potential, galerkin);
    else
      std::cout << input.controlVariateName() << " is not a valid control variate !" << std::endl;
    return 0;
  }


  ControlVariate::ControlVariate(Input const& input, const string& outPath, Potential& potential, int nbOfFunctions):
    Observable(input, outPath),
    dimension_(input.dimension()),
    //decorrelationNbOfSteps_(input.decorrelationNbOfSteps()),
    //timeStep_(input.timeStep()),
    nbOfFunctions_(nbOfFunctions),
    nbOfFunctionPairs_(pow(nbOfFunctions_, 2)),
    //printPeriodNbOfSteps_(input.printPeriodNbOfSteps()),
    //nbOfAutocoPts_(input.nbOfAutocoPts()),
    //autocoStats_(decorrelationNbOfSteps(), timeStep(), nbOfAutocoPts()),
    autocoStatsBetter_(decorrelationNbOfSteps(), timeStep(), nbOfAutocoPts()),
    statsGeneratorOnBasis_(nbOfFunctions_),
    statsB1_(nbOfFunctions_),
    statsB2_(decorrelationNbOfSteps(), timeStep(), nbOfAutocoPts(), nbOfFunctions_),
    statsD_(nbOfFunctions_, nbOfFunctions_),
    lastA_(nbOfFunctions_),
    //historyObservable_(input.nbOfSteps()),
    //historyGeneratorOnBasis_(input.nbOfSteps(), nbOfFunctions_),
    potential_(&potential),
    valueBasis_(nbOfFunctions_),
    valueGeneratorOnBasis_(nbOfFunctions_)
  {}

  int ControlVariate::nbOfFunctions() const
  {
    return nbOfFunctions_;
  }

  int ControlVariate::nbOfFunctionPairs() const
  {
    return nbOfFunctionPairs_;
  }

  /*bool ControlVariate::doOutput(long int iOfStep) const
  {
    return (printPeriodNbOfSteps_ > 0 && iOfStep % printPeriodNbOfSteps_ == 0);
  }*/

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
    return potential_->gradient(position);
  }

  double ControlVariate::potentialLaplacian(Vector<double> const& position) const
  {
    return potential_->laplacian(position);
  }

  /*const int& ControlVariate::decorrelationNbOfSteps() const
  {
    return decorrelationNbOfSteps_;
  }
  
  const double& ControlVariate::timeStep() const
  {
    return timeStep_;
  }

  double ControlVariate::decorrelationTime() const
  {
    return decorrelationNbOfSteps() * timeStep();
  }

  int const& ControlVariate::nbOfAutocoPts() const
  {
    return nbOfAutocoPts_;
  }*/

  int ControlVariate::nbOfFourier() const
  {return 0;}

  int ControlVariate::nbOfHermite() const
  {return 0;}




  /*double ControlVariate::lastObservable() const
  {
    return autocoStats_.lastValue();
  }


  double ControlVariate::meanObservable() const
  {
    return autocoStats_.mean();
  }
  
  double ControlVariate::varObservable() const
  {
    return autocoStats_.variance();
  }

  double ControlVariate::stdDevObservable() const
  {
    return autocoStats_.stdDev();
  }
  
  double ControlVariate::varOfVarObservable() const
  {
    return autocoStats_.varOfVar();
  }
  
  double ControlVariate::stdDevOfVarObservable() const
  {
    return autocoStats_.stdDevOfVar();
  }*/

  double ControlVariate::lastValueBetter() const
  {
    return autocoStatsBetter_.lastValue();
  }

  double ControlVariate::meanBetter() const
  {
    return autocoStatsBetter_.mean();
  }
  
  double ControlVariate::varBetter() const
  {
    return autocoStatsBetter_.variance();
  }

  double ControlVariate::stdDevBetter() const
  {
    return autocoStatsBetter_.stdDev();
  }
  
  double ControlVariate::varOfVarBetter() const
  {
    return autocoStatsBetter_.varOfVar();
  }

  Vector<double> ControlVariate::correlationB2() const
  {
    return statsB2_.integratedCorrelationVec();
  }

  double ControlVariate::correlationB2(int iOfFunction) const
  {
    return statsB2_.integratedCorrelation(iOfFunction);
  }

  void ControlVariate::appendToObservable(double observable, long int iOfStep)
  {
    //historyObservable_(iOfStep) = observable;
    autocoStats_.append(observable, iOfStep);
  }

  /*double ControlVariate::correlationAtSpan(int iOfSpan) const
  {
    return autocoStats_.correlationAtSpan(iOfSpan);
  }
  
  double ControlVariate::unbiasedCorrelationAtSpan(int iOfSpan) const
  {
    //return autocoStats_(iOfSpan) - pow(meanObservable(), 2);
    return autocoStats_.unbiasedCorrelationAtSpan(iOfSpan);
  }
  
  double ControlVariate::stdDevCorrelationAtSpan(int iOfSpan) const
  {
    return autocoStats_.stdDevCorrelationAtSpan(iOfSpan);
  }
  
  double ControlVariate::varCorrelationAtSpan(int iOfSpan) const
  {
    return autocoStats_.varCorrelationAtSpan(iOfSpan);
  }*/
  
  /*double ControlVariate::stdErrorCorrelationAtSpan(int iOfSpan) const
  {
    return autocoStats_.stdErrorCorrelationAtSpan(iOfSpan);
  }*/

  void ControlVariate::appendToB1(double observable)
  {
    for (int iOfFunction = 0; iOfFunction < nbOfFunctions_; iOfFunction++)
      statsB1_.append((observable - autocoStats_.mean()) * valueBasis_(iOfFunction), iOfFunction);
  }

  void ControlVariate::appendToB2(double observable, long int iOfStep)
  {
    for (int iOfFunction = 0; iOfFunction < nbOfFunctions_; iOfFunction++)
      statsB2_.append((observable - autocoStats_.mean()), iOfStep, iOfFunction, valueGeneratorOnBasis_(iOfFunction));
  }

  void ControlVariate::appendToD()
  {
    //Eigen::Matrix<CorrelationStats<double>, Eigen::Dynamic, Eigen::Dynamic> A(nbOfFunctions_, nbOfFunctions_, CorrelationStats<double>(decorrelationNbOfSteps(), decorrelationTime()));

    for (int iOfFunction = 0; iOfFunction < nbOfFunctions_; iOfFunction++)
      for (int iOfFunction2 = 0; iOfFunction2 <= iOfFunction; iOfFunction2++)
      {
        double valueSym = (- valueBasis_(iOfFunction) * valueGeneratorOnBasis_(iOfFunction2)
                           - valueBasis_(iOfFunction2) * valueGeneratorOnBasis_(iOfFunction)) / 2.;
        statsD_.append(valueSym, iOfFunction, iOfFunction2);
        if (iOfFunction != iOfFunction2)
          statsD_.append(valueSym, iOfFunction2, iOfFunction);
      }
  }



  void ControlVariate::appendToBetter(double observable, long int iOfStep)
  {
    //double betterObservableTerm = dot((statsB_.meanMat() - statsB2_.integratedCorrelationMat()), statsD_.meanMat().llt().solve(generatorOnBasisFunction));
    /*DenseMatrix<double> Dinv = statsD_.meanMat().llt().solve(generatorOnBasisFunction);
    Vector<double> B = statsB_.meanMat() - statsB2_.integratedCorrelationMat();
    std::cout << B.transpose() * Dinv << endl;*/
    //cout << (statsB_.meanMat() - statsB2_.integratedCorrelationMat()).transpose().size() << "  " << statsD_.meanMat().inverse().size() << endl;
    //lastA_ = (statsB_.meanMat() - statsB2_.integratedCorrelationMat()).transpose() * statsD_.meanMat().llt().solve(generatorOnBasisFunction);
    //cout << lastA_.size() << "   " << generatorOnBasisFunction.size() << endl;

    if (statsD_.meanMat().determinant() != 0)
    {
      //lastA_ = - .5 * statsD_.meanMat().llt().solve(meanB());
      //lastA_ = Vector<double>(1,1);
      lastA_.fill(1);
      autocoStatsBetter_.append(observable - dot(lastA_, valueGeneratorOnBasis_), iOfStep);
    }
    else
    {
      lastA_.fill(0);
      autocoStatsBetter_.append(observable, iOfStep);
    }

    for (int iOfFunction = 0; iOfFunction < nbOfFunctions_; iOfFunction++)
    {
      //historyGeneratorOnBasis_(iOfStep, iOfFunction) = generatorOnBasisFunction(iOfFunction);
      statsGeneratorOnBasis_.append(valueGeneratorOnBasis_(iOfFunction), iOfFunction);
    }
  }



  Vector<double> ControlVariate::lastB1() const
  {
    /*Vector<double> result(nbOfFunctions_);
    for (int iOfFunction =0; iOfFunction < nbOfFunctions_; iOfFunction++)
      result(iOfFunction) = statsB_.lastValue(iOfFunction);
    return result;*/
    return statsB1_.lastValueVec();
  }


  Vector<double> ControlVariate::meanB1() const
  {
    return statsB1_.meanVec();
  }

  Vector<double> ControlVariate::meanB() const
  {
    return statsB1_.meanVec() - statsB2_.integratedCorrelationVec();
  }

  DenseMatrix<double> ControlVariate::lastD() const
  {
    return statsD_.lastValueMat();
  }

  DenseMatrix<double> ControlVariate::meanD() const
  {
    return statsD_.meanMat();
  }

  Vector<double> ControlVariate::lastGeneratorOnBasis() const
  {
    return statsGeneratorOnBasis_.lastValueVec();
  }

  Vector<double> ControlVariate::meanGeneratorOnBasis() const
  {
    return statsGeneratorOnBasis_.meanVec();
  }

  double ControlVariate::autocorrelationB2(int iOfSpan, int iOfFunction) const
  {
    return statsB2_.correlationAtSpan(iOfSpan, iOfFunction);
  }

  Vector<double> ControlVariate::lastA() const
  {
    return lastA_;

  }
  double ControlVariate::lastA(int iOfFunction) const
  {
    return lastA_(iOfFunction);
  }

  ///
  ///Feed a new value to the Observable
  /// /!\ an update*** method must be called beforehand
  void ControlVariate::append(double observable, long int iOfStep)
  {
    appendToB1(observable);
    appendToB2(observable, iOfStep);
    appendToD();
    appendToObservable(observable, iOfStep);
    appendToBetter(observable, iOfStep);
    //cout << "end ControlVariate::update" << endl;
  }
  
  void ControlVariate::updateHamiltonian(vector<Particle> const& configuration)
  {
    computeValueBasis(configuration);
    computeGeneratorHamiltonian(configuration);
  }
  
  void ControlVariate::updateOverdamped(vector<Particle> const& configuration, double beta)
  {
    computeValueBasis(configuration);
    computeGeneratorOverdamped(configuration, beta);
  }
  
  void ControlVariate::updateLangevin(vector<Particle> const& configuration, double beta, double gamma)
  {
    computeValueBasis(configuration);
    computeGeneratorLangevin(configuration, beta, gamma);
  }
  
  void ControlVariate::updateBoundaryLangevin(vector<Particle> const& configuration, double betaLeft, double betaRight, double gamma)
  {
    computeValueBasis(configuration);
    computeGeneratorBoundaryLangevin(configuration, betaLeft, betaRight, gamma);
  }
  
  ///
  ///Evaluates each function of the basis at "configuration" and updates the attribute
  void ControlVariate::computeValueBasis(vector<Particle> const& configuration)
  {
    for (int iOfFunction = 0; iOfFunction < nbOfFunctions_; iOfFunction++)
      valueBasis_(iOfFunction) = basisFunction(configuration, iOfFunction);
  }

  ///
  ///Applies the generator of this dynamics to the basis functions of the CV
  void ControlVariate::computeGeneratorHamiltonian(vector<Particle> const& configuration)
  {
    int nbOfParticles = (int)configuration.size();
    //Vector<double> result = Vector<double>::Zero(nbOfFunctions());
    for (int iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
      for (int iOfParticle = 0; iOfParticle < nbOfParticles; iOfParticle++)
        valueGeneratorOnBasis_(iOfFunction) += dot( configuration[iOfParticle].momentum() , gradientQ(configuration, iOfParticle, iOfFunction))
                               + dot( configuration[iOfParticle].force() , gradientP(configuration, iOfParticle, iOfFunction));
  }

  ///Applies the generator of this dynamics to the basis functions of the CV
  ///Evaluate at the current state of "conifguration"
  void ControlVariate::computeGeneratorOverdamped(vector<Particle> const& configuration, double beta)
  {
    int nbOfParticles = (int)configuration.size();
    //Vector<double> result = Vector<double>::Zero(nbOfFunctions());
    for (int iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
      for (int iOfParticle = 0; iOfParticle < nbOfParticles; iOfParticle++)
        valueGeneratorOnBasis_(iOfFunction) += laplacianQ(configuration, iOfParticle, iOfFunction) / beta
                               + dot(configuration[iOfParticle].force(), gradientQ(configuration, iOfParticle, iOfFunction));
  }

  ///Applies the generator of this dynamics to the basis functions of the CV
  ///Evaluate at the current state of "conifguration"
  void ControlVariate::computeGeneratorLangevin(vector<Particle> const& configuration, double beta, double gamma)
  {
    int nbOfParticles = (int)configuration.size();
    //cout << "generatorOn(const Langevin& dyna, S const& syst, const ControlVariate& controlVariate)" << endl;
    //Vector<double> result = Vector<double>::Zero(nbOfFunctions());
    for (int iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
      for (int iOfParticle = 0; iOfParticle < nbOfParticles; iOfParticle++)
      {
        valueGeneratorOnBasis_(iOfFunction) += dot(configuration[iOfParticle].momentum() , gradientQ(configuration, iOfParticle, iOfFunction))
                               + dot(configuration[iOfParticle].force() , gradientP(configuration, iOfParticle, iOfFunction))
                               + gamma * (- dot(configuration[iOfParticle].momentum() , gradientP(configuration, iOfParticle, iOfFunction))
                                          + laplacianP(configuration, iOfParticle, iOfFunction) / beta );
      }
  }

  ///Applies the generator of this dynamics to the basis functions of the CV
  ///Evaluate at the current state of "conifguration"
  void ControlVariate::computeGeneratorBoundaryLangevin(vector<Particle> const& configuration, double betaLeft, double betaRight, double gamma)
  {
    int nbOfParticles = (int)configuration.size();
    //Vector<double> result = Vector<double>::Zero(nbOfFunctions());
    for (int iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
    {
      for (int iOfParticle = 0; iOfParticle < nbOfParticles; iOfParticle++)
        valueGeneratorOnBasis_(iOfFunction) += dot(configuration[iOfParticle].momentum(), gradientQ(configuration, iOfParticle, iOfFunction))
                               + dot(configuration[iOfParticle].force(), gradientP(configuration, iOfParticle, iOfFunction));
      //if(false)
      valueGeneratorOnBasis_(iOfFunction) += gamma * (- dot(configuration[0].momentum(), gradientP(configuration, 0, iOfFunction))
                                      + laplacianP(configuration, 0, iOfFunction) / betaLeft
                                      - dot(configuration[nbOfParticles - 1].momentum(), gradientP(configuration, nbOfParticles - 1, iOfFunction))
                                      + laplacianP(configuration, nbOfParticles - 1, iOfFunction) / betaRight);
    }
  }


  void ControlVariate::display(std::ofstream& out, double time) const
  {
    out << time;
    for (int iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
      out << " " << meanB()(iOfFunction) * lastA(iOfFunction);
    for (int iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
      out << " " << lastA(iOfFunction);
    /*for (int iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
      out << " " << lastB()(iOfFunction);*/
    for (int iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
      out << " " << meanB()(iOfFunction);
    /*for (int iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
      for (int iOfFunction2 = 0; iOfFunction2 < nbOfFunctions(); iOfFunction2++)
    out << " " << lastD()(iOfFunction, iOfFunction2);*/
    for (int iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
      for (int iOfFunction2 = 0; iOfFunction2 < nbOfFunctions(); iOfFunction2++)
        out << " " << meanD()(iOfFunction, iOfFunction2);  //8-11

    out << " " << lastValue()
        << " " << mean()   //13
        << " " << variance()
        << " " << varOfVar()
        << " " << lastValueBetter()
        << " " << meanBetter()
        << " " << varOfVarBetter();
    for (int iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
      out << " " << lastGeneratorOnBasis()(iOfFunction);
    for (int iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
      out << " " << meanGeneratorOnBasis()(iOfFunction);

    for (int iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
      out << " " << meanB1()(iOfFunction);

    for (int iOfFunction = 0; iOfFunction < nbOfFunctions(); iOfFunction++)
      out << " " << -correlationB2()(iOfFunction);  //22-23

    out << endl;
  }



 

  SinusControlVariate::SinusControlVariate(Input const& input, const string& outPath, Potential& potential):
    ControlVariate(input, outPath, potential, 1)
  {}

  double SinusControlVariate::basisFunction(vector<Particle> const& configuration, int /*iOfFunction*/) const
  {
    double q = configuration[0].position(0);
    return sin(2 * M_PI * q);
  }

  Vector<double> SinusControlVariate::gradientQ(vector<Particle> const& configuration, int /*iOfParticle*/, int /*iOfFunction*/) const
  {
    double q = configuration[0].position(0);
    Vector<double> grad(dimension_);
    grad(0) = 2 * M_PI * cos(2 * M_PI * q);
    return grad;
  }

  double SinusControlVariate::laplacianQ(vector<Particle> const& configuration, int /*iOfParticle*/, int /*iOfFunction*/) const
  {
    double q = configuration[0].position(0);
    return - pow (2 * M_PI, 2) * sin(2 * M_PI * q);
  }


  double SinusControlVariate::laplacianP(vector<Particle> const& /*configuration*/, int /*iOfParticle*/, int /*iOfFunction*/) const
  {
    return 0;
  }

  Vector<double> SinusControlVariate::gradientP(vector<Particle> const& /*configuration*/, int /*iOfParticle*/, int /*iOfFunction*/) const
  {
    return Vector<double>(dimension_);
  }




  CosControlVariate::CosControlVariate(Input const& input, const string& outPath, Potential& potential):
    ControlVariate(input, outPath, potential, 1)
  {}

  double CosControlVariate::basisFunction(vector<Particle> const& configuration, int /*iOfFunction*/) const
  {
    double q = configuration[0].position(0);
    return cos(2 * M_PI * q);
  }

  Vector<double> CosControlVariate::gradientQ(vector<Particle> const& configuration, int /*iOfParticle*/, int /*iOfFunction*/) const
  {
    double q = configuration[0].position(0);
    Vector<double> grad(dimension_);
    grad(0) = - 2 * M_PI * sin(2 * M_PI * q);
    return grad;
  }

  double CosControlVariate::laplacianQ(vector<Particle> const& configuration, int /*iOfParticle*/, int /*iOfFunction*/) const
  {
    double q = configuration[0].position(0);
    return - pow (2 * M_PI, 2) * cos(2 * M_PI * q);
  }

  Vector<double> CosControlVariate::gradientP(vector<Particle> const& /*configuration*/, int /*iOfParticle*/, int /*iOfFunction*/) const
  {
    return Vector<double>(dimension_);
  }

  double CosControlVariate::laplacianP(vector<Particle> const& /*configuration*/, int /*iOfParticle*/, int /*iOfFunction*/) const
  {
    return 0;
  }








  SinExpControlVariate::SinExpControlVariate(Input const& input, const string& outPath, Potential& potential):
    ControlVariate(input, outPath, potential, 1)
  {}

  double SinExpControlVariate::basisFunction(vector<Particle> const& configuration, int /*iOfFunction*/) const
  {
    Vector<double> q = configuration[0].position();
    double q0 = q(0);
    return sin(2 * M_PI * q0) * exp(potential(q) / 2);
  }

  Vector<double> SinExpControlVariate::gradientQ(vector<Particle> const& configuration, int /*iOfParticle*/, int /*iOfFunction*/) const
  {
    Vector<double> q = configuration[0].position();
    double q0 = q(0);
    Vector<double> grad(dimension_);
    grad(0) = (2 * M_PI * cos(2 * M_PI * q0)
               + potentialDerivative(q)(0) / 2 * sin(2 * M_PI * q0))
              * exp(potential(q) / 2);
    return grad;
  }

  double SinExpControlVariate::laplacianQ(vector<Particle> const& configuration, int /*iOfParticle*/, int /*iOfFunction*/) const
  {
    Vector<double> q = configuration[0].position();
    double q0 = q(0);

    return (-pow(2 * M_PI, 2) * sin(2 * M_PI * q0)
            + 2 * M_PI * potentialDerivative(q)(0) *  cos(2 * M_PI * q0)
            + potentialLaplacian(q) / 2 * sin(2 * M_PI * q0)
            + pow(potentialDerivative(q)(0) / 2, 2) * sin(2 * M_PI * q0))
           * exp(potential(q) / 2);
  }

  Vector<double> SinExpControlVariate::gradientP(vector<Particle> const& /*configuration*/, int /*iOfParticle*/, int /*iOfFunction*/) const
  {
    return Vector<double>(dimension_);
  }

  double SinExpControlVariate::laplacianP(vector<Particle> const& /*configuration*/, int /*iOfParticle*/, int /*iOfFunction*/) const
  {
    return 0;
  }




  CosExpControlVariate::CosExpControlVariate(Input const& input, const string& outPath, Potential& potential):
    ControlVariate(input, outPath, potential, 1)
  {}

  double CosExpControlVariate::basisFunction(vector<Particle> const& configuration, int /*iOfFunction*/) const
  {
    Vector<double> q = configuration[0].position();
    double q0 = q(0);
    return cos(2 * M_PI * q0) * exp(potential(q) / 2);
  }

  Vector<double> CosExpControlVariate::gradientQ(vector<Particle> const& configuration, int /*iOfParticle*/, int /*iOfFunction*/) const
  {
    Vector<double> q = configuration[0].position();
    double q0 = q(0);
    Vector<double> grad(dimension_);
    grad(0) = (- 2 * M_PI * sin(2 * M_PI * q0)
               + potentialDerivative(q)(0) / 2 * cos(2 * M_PI * q0))
              * exp(potential(q) / 2);
    return grad;
  }

  double CosExpControlVariate::laplacianQ(vector<Particle> const& configuration, int /*iOfParticle*/, int /*iOfFunction*/) const
  {
    Vector<double> q = configuration[0].position();
    double q0 = q(0);

    return (-pow(2 * M_PI, 2) * cos(2 * M_PI * q0)
            - 2 * M_PI * potentialDerivative(q)(0) *  sin(2 * M_PI * q0)
            + potentialLaplacian(q) / 2 * cos(2 * M_PI * q0)
            + pow(potentialDerivative(q)(0) / 2, 2) * cos(2 * M_PI * q0))
           * exp(potential(q) / 2);
  }

  Vector<double> CosExpControlVariate::gradientP(vector<Particle> const& /*configuration*/, int /*iOfParticle*/, int /*iOfFunction*/) const
  {
    return Vector<double>(dimension_);
  }

  double CosExpControlVariate::laplacianP(vector<Particle> const& /*configuration*/, int /*iOfParticle*/, int /*iOfFunction*/) const
  {
    return 0;
  }




  LangevinControlVariate::LangevinControlVariate(Input const& input, const string& outPath, Potential& potential):
    ControlVariate(input, outPath, potential, 1)
  {}

  double LangevinControlVariate::basisFunction(vector<Particle> const& configuration, int /*iOfFunction*/) const
  {
    //return pow(configuration[0].momentum(0), 2) * configuration[0].position(0);
    //return configuration[0].momentum(0);
    return configuration[0].kineticEnergy();
  }

  Vector<double> LangevinControlVariate::gradientQ(vector<Particle> const& /*configuration*/, int /*iOfParticle*/, int /*iOfFunction*/) const
  {
    Vector<double> grad(dimension_);
    //grad(0) = pow(configuration[0].momentum(0), 2);
    grad(0) = 0;
    return grad;
  }

  double LangevinControlVariate::laplacianQ(vector<Particle> const& /*configuration*/, int /*iOfParticle*/, int /*iOfFunction*/) const
  {
    throw std::runtime_error("LangevinControlVariate::laplacianQ should not be called");
    return 0;
  }

  Vector<double> LangevinControlVariate::gradientP(vector<Particle> const& configuration, int /*iOfParticle*/, int /*iOfFunction*/) const
  {
    Vector<double> grad(dimension_);
    grad(0) = configuration[0].momentum(0);
    //grad(0) = 1;
    return grad;
  }

  double LangevinControlVariate::laplacianP(vector<Particle> const& /*configuration*/, int /*iOfParticle*/, int /*iOfFunction*/) const
  {
    return 1;
    //return 2 * configuration[0].position(0);
  }



  SumEnergyControlVariate::SumEnergyControlVariate(Input const& input, const string& outPath, Potential& potential):
    ControlVariate(input, outPath, potential, 1), i0_(0)
  {}

  double SumEnergyControlVariate::basisFunction(vector<Particle> const& configuration, int /*iOfFunction*/) const
  {
    double result = 0;
    for (int iOfParticle = 1; iOfParticle < (int)configuration.size(); iOfParticle++)
      result += configuration[iOfParticle].kineticEnergy() * (iOfParticle - i0_)
                + (iOfParticle - i0_ - .5) * configuration[iOfParticle].potentialEnergy();
    //cout << result << endl;
    return result;
  }

  Vector<double> SumEnergyControlVariate::gradientQ(vector<Particle> const& configuration, int iOfParticle, int /*iOfFunction*/) const
  {
    return (iOfParticle - i0_ - .5) * configuration[iOfParticle].energyGrad();
  }

  double SumEnergyControlVariate::laplacianQ(vector<Particle> const& /*configuration*/, int /*iOfParticle*/, int /*iOfFunction*/) const
  {
    throw std::runtime_error("SumEnergyControlVariate::laplacianQ should not be called");
    return 0;
  }

  Vector<double> SumEnergyControlVariate::gradientP(vector<Particle> const& configuration, int iOfParticle, int /*iOfFunction*/) const
  {
    return (iOfParticle - i0_) * configuration[iOfParticle].momentum();
  }

  double SumEnergyControlVariate::laplacianP(vector<Particle> const& /*configuration*/, int iOfParticle, int /*iOfFunction*/) const
  {
    return iOfParticle - i0_;
  }




  EnergyControlVariate::EnergyControlVariate(Input const& input, const string& outPath, Potential& potential):
    ControlVariate(input, outPath, potential, 1), i0_(0)
  {}

  double EnergyControlVariate::basisFunction(vector<Particle> const& configuration, int /*iOfFunction*/) const
  {
    double result = 0;
    for (int iOfParticle = 0; iOfParticle < (int)configuration.size(); iOfParticle++)
      result += configuration[iOfParticle].kineticEnergy()
                + configuration[iOfParticle].potentialEnergy();
    //cout << result << endl;
    return result;
  }

  Vector<double> EnergyControlVariate::gradientQ(vector<Particle> const& configuration, int iOfParticle, int /*iOfFunction*/) const
  {
    return configuration[iOfParticle].energyGrad();
  }

  double EnergyControlVariate::laplacianQ(vector<Particle> const& /*configuration*/, int /*iOfParticle*/, int /*iOfFunction*/) const
  {
    throw std::runtime_error("EnergyControlVariate::laplacianQ should not be called");
  }

  Vector<double> EnergyControlVariate::gradientP(vector<Particle> const& configuration, int iOfParticle, int /*iOfFunction*/) const
  {
    return configuration[iOfParticle].momentum();
  }

  double EnergyControlVariate::laplacianP(vector<Particle> const& /*configuration*/, int /*iOfParticle*/, int /*iOfFunction*/) const
  {
    return 1;
  }



  LocalControlVariate::LocalControlVariate(Input const& input, const string& outPath, Potential& potential):
    ControlVariate(input, outPath, potential, 1)
  {}

  double LocalControlVariate::basisFunction(vector<Particle> const& configuration, int /*iOfFunction*/) const
  {
    //return configuration[0].kineticEnergy() + configuration[0].potentialEnergy();
    return configuration[0].energy();
  }

  Vector<double> LocalControlVariate::gradientQ(vector<Particle> const& configuration, int iOfParticle, int /*iOfFunction*/) const
  {
    Vector<double> result = Vector<double>(dimension_);

    if (iOfParticle == 0)
      result(0) = sin(configuration[0].position(0));

    return result;
  }

  double LocalControlVariate::laplacianQ(vector<Particle> const& /*configuration*/, int /*iOfParticle*/, int /*iOfFunction*/) const
  {
    throw std::runtime_error("LocalControlVariate::laplacianQ should not be called");
  }

  Vector<double> LocalControlVariate::gradientP(vector<Particle> const& configuration, int iOfParticle, int /*iOfFunction*/) const
  {
    Vector<double> result = Vector<double>(dimension_);

    if (iOfParticle == 0)
      result(0) = configuration[0].momentum(0);

    return result;
  }

  double LocalControlVariate::laplacianP(vector<Particle> const& /*configuration*/, int iOfParticle, int /*iOfFunction*/) const
  {
    if (iOfParticle == 0)
      return 1;
    else
      return 0;
  }



  KineticControlVariate::KineticControlVariate(Input const& input, const string& outPath, Potential& potential):
    ControlVariate(input, outPath, potential, 1)
  {}

  double KineticControlVariate::basisFunction(vector<Particle> const& configuration, int /*iOfFunction*/) const
  {
    return 2 * configuration[0].momentum(0) * (sin(configuration[1].position(0) - configuration[0].position(0)) - sin(configuration[0].position(0)))
           + 2 * configuration[0].kineticEnergy();
  }

  Vector<double> KineticControlVariate::gradientQ(vector<Particle> const& configuration, int iOfParticle, int /*iOfFunction*/) const
  {
    Vector<double> result(dimension_);

    if (iOfParticle == 0)
      result(0) = -2 * configuration[0].momentum(0) * (cos(configuration[1].position(0) - configuration[0].position(0)) + cos(configuration[0].position(0)));
    else if (iOfParticle == 1)
      result(0) = 2 * configuration[0].momentum(0) * cos(configuration[1].position(0) - configuration[0].position(0));

    return result;
  }

  double KineticControlVariate::laplacianQ(vector<Particle> const& /*configuration*/, int /*iOfParticle*/, int /*iOfFunction*/) const
  {
    throw std::runtime_error("KineticControlVariate::laplacianQ should not be called");
  }

  Vector<double> KineticControlVariate::gradientP(vector<Particle> const& configuration, int iOfParticle, int /*iOfFunction*/) const
  {
    Vector<double> result(dimension_);

    if (iOfParticle == 0)
      result(0) = 2 * (sin(configuration[1].position(0) - configuration[0].position(0)) - sin(configuration[0].position(0)))
                  + 2 * configuration[0].momentum(0);

    return result;
  }

  double KineticControlVariate::laplacianP(vector<Particle> const& /*configuration*/, int iOfParticle, int /*iOfFunction*/) const
  {
    if (iOfParticle == 0)
      return 2;
    else
      return 0;
  }



  //#####TwoControlVariate#####

  TwoControlVariate::TwoControlVariate(Input const& input, const string& outPath, Potential& potential):
    ControlVariate(input, outPath, potential, 2)
  {}

  double TwoControlVariate::basisFunction(vector<Particle> const& configuration, int iOfFunction) const
  {
    if (iOfFunction == 0)
      return (configuration[0].momentum(0) - configuration[configuration.size() - 1].momentum(0))
             * (configuration[0].position(0) + configuration[configuration.size() - 1].position(0));
    else if (iOfFunction == 1)
      return (configuration[0].momentum(0) + configuration[configuration.size() - 1].momentum(0))
             * (configuration[0].position(0) - configuration[configuration.size() - 1].position(0));
    else
      throw std::invalid_argument("iOfFunction must 0 or 1");
  }

  Vector<double> TwoControlVariate::gradientQ(vector<Particle> const& configuration, int iOfParticle, int iOfFunction) const
  {
    Vector<double> result = Vector<double>(dimension_);

    if (iOfFunction == 0)
      result(0) = - configuration[0].momentum(0) - configuration[configuration.size() - 1].momentum(0);
    else if (iOfFunction == 1)
    {
      if (iOfParticle == 0)
        result(0) = configuration[0].momentum(0) + configuration[configuration.size() - 1].momentum(0);
      else if (iOfParticle == 1)
        result(0) = - configuration[0].momentum(0) - configuration[configuration.size() - 1].momentum(0);
    }
    return result;
  }

  double TwoControlVariate::laplacianQ(vector<Particle> const& /*configuration*/, int /*iOfParticle*/, int /*iOfFunction*/) const
  {
    throw std::runtime_error("TwoControlVariate::laplacianQ should not be called");
  }

  Vector<double> TwoControlVariate::gradientP(vector<Particle> const& configuration, int iOfParticle, int iOfFunction) const
  {
    Vector<double> result(dimension_);

    if (iOfFunction == 0)
    {
      if (iOfParticle == 0)
        result(0) = configuration[0].position(0) + configuration[configuration.size() - 1].position(0);
      else if (iOfParticle == 1)
        result(0) = - configuration[0].position(0) - configuration[configuration.size() - 1].position(0);
    }
    else if (iOfFunction == 1)
      result(0) = - configuration[0].position(0) - configuration[configuration.size() - 1].position(0);


    return result;
  }

  double TwoControlVariate::laplacianP(vector<Particle> const& /*configuration*/, int /*iOfParticle*/, int /*iOfFunction*/) const
  {
    return 0;
  }


  //#####BasisControlVariate#####

  BasisControlVariate::BasisControlVariate(Input const& input, const string& outPath, Potential& potential, Galerkin* galerkin):
    ControlVariate(input, outPath, potential, 1),
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
      std::string coeffsPath = input.outputFolderName() + input.controlVariateCoeffsPath();
      ifstream file(coeffsPath);
      coeffsVec_ = SparseMatrix<double>(coeffsPath, getNbOfLines(file));
    }
    //cout << "coeffsVec_ : " << endl;
    //cout << coeffsVec_ << endl;
  }

  double BasisControlVariate::basisFunction(vector<Particle> const& configuration, int iOfFunction) const
  {
    assert(iOfFunction == 0);
    double result = 0;

    for (SMat::iterator it(coeffsVec_, 0); it; ++it)
    {
      int iOfCoeff = it.row();
      double valOfCoeff = it.value();
      result += valOfCoeff * basis_->value(configuration, iOfCoeff);
    }

    return result;
  }

  Vector<double> BasisControlVariate::gradientQ(vector<Particle> const& configuration, int iOfParticle, int iOfFunction) const
  {
    //cout << "BasisControlVariate::gradientQ" << endl;
    Vector<double> result(1, 0);
    assert(iOfFunction == 0);
    //cout << "gradientQ" << endl;
    for (SMat::iterator it(coeffsVec_, 0); it; ++it)
    {
      int iOfCoeff = it.row();
      double valOfCoeff = it.value();
      //cout << "+ " << valOfCoeff << " X ";
      result += valOfCoeff * basis_->gradientQ(configuration, iOfParticle, iOfCoeff);
    }
    //cout << endl;
    //cout << "end BasisControlVariate::gradientQ" << endl;
    //cout << "gradientQ = " << result << endl;
    return result;
  }

  double BasisControlVariate::laplacianQ(vector<Particle> const& configuration, int iOfParticle, int iOfFunction) const
  {
    double result = 0;
    assert(iOfFunction == 0);

    for (SMat::iterator it(coeffsVec_, 0); it; ++it)
    {
      int iOfCoeff = it.row();
      double valOfCoeff = it.value();
      result += valOfCoeff * basis_->laplacianQ(configuration, iOfParticle, iOfCoeff);
    }
    //cout << "laplacianQ = " << result << endl;
    return result;
  }

  Vector<double> BasisControlVariate::gradientP(vector<Particle> const& configuration, int iOfParticle, int iOfFunction) const
  {
    Vector<double> result(1, 0);
    assert(iOfFunction == 0);
    //cout << "gradientP" << endl;
    for (SMat::iterator it(coeffsVec_, 0); it; ++it)
    {
      int iOfCoeff = it.row();
      double valOfCoeff = it.value();
      //cout << "+ " << valOfCoeff << " X ";
      result += valOfCoeff * basis_->gradientP(configuration, iOfParticle, iOfCoeff);
    }
    //cout << endl;
    //cout << "gradientP = " << result << endl;
    return result;
  }

  double BasisControlVariate::laplacianP(vector<Particle> const& configuration, int iOfParticle, int iOfFunction) const
  {
    double result = 0;
    assert(iOfFunction == 0);

    for (SMat::iterator it(coeffsVec_, 0); it; ++it)
    {
      int iOfCoeff = it.row();
      double valOfCoeff = it.value();
      result += valOfCoeff * basis_->laplacianP(configuration, iOfParticle, iOfCoeff);
    }
    //cout << "laplacianP = " << result << endl;
    return result;
  }



  ExpFourierHermiteControlVariate::ExpFourierHermiteControlVariate(Input const& input, const string& outPath, Potential& potential, Galerkin* galerkin):
    BasisControlVariate(input, outPath, potential, galerkin)
  {
    basis_ = new ExpFourierHermiteBasis(input, potential);
    nbQ_ = 100;  // rafinement du maillage de l'output map
    nbP_ = 100;  // rafinement du maillage de l'output map
    deltaQ_ = 2 * M_PI / nbQ_;
    pMax_ = 4;
    deltaP_ = 2 * pMax_ / nbP_;
  }

  int ExpFourierHermiteControlVariate::nbOfFourier() const
  {return basis_->nbOfElts(0);}

  int ExpFourierHermiteControlVariate::nbOfHermite() const
  {return basis_->nbOfElts(1);}

  void ExpFourierHermiteControlVariate::displayMap(ofstream& out) const
  {
    cout << "displayMap" << endl;
    vector<Particle> conf(1, Particle(0, 0, 0));
    //conf[0] = Particle(0,0,0);
    out << nbP_ + 1 << " ";
    for (double p = -pMax_; p < pMax_; p += deltaP_)
      out << p << " ";
    out << endl;
    for (double q = -M_PI; q < M_PI; q += deltaQ_)
    {
      out << q << " ";
      for (double p = -pMax_; p < pMax_; p += deltaP_)
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
    vector<Particle> conf(1, Particle(0, 0, 0));
    //conf[0] = Particle(0,0,0);
    out << nbP_ + 1 << " ";
    for (double p = -pMax_; p < pMax_; p += deltaP_)
      out << p << " ";
    out << endl;
    for (double q = -M_PI; q < M_PI; q += deltaQ_)
    {
      out << q << " ";
      for (double p = -pMax_; p < pMax_; p += deltaP_)
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
    vector<Particle> conf(1, Particle(0, 0, 0));
    //nf[0] = Particle(0,0,0);
    out << nbP_ + 1 << " ";
    for (double p = -pMax_; p < pMax_; p += deltaP_)
      out << p << " ";
    out << endl;
    for (double q = -M_PI; q < M_PI; q += deltaQ_)
    {
      out << q << " ";
      for (double p = -pMax_; p < pMax_; p += deltaP_)
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










