#include "controlVariate.hpp"

using std::cout; 
using std::endl; 
using std::sin;
using std::cos;
using std::exp;

namespace simol 
{
  
  ControlVariate* createControlVariate(Input const& input, Potential& potential, size_t indexOfReplica)
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
    else
      std::cout << input.controlVariateName() << " is not a valid control variate !" << std::endl;
    return 0;
  }

  
  ControlVariate::ControlVariate(Input const& input, Potential& potential):
    statsB_(input.decorrelationNumberOfIterations()),
    statsD_(input.decorrelationNumberOfIterations()),
    statsObservable_(input.decorrelationNumberOfIterations()),
    statsBetterObservable_(input.decorrelationNumberOfIterations()),
    statsGeneratorOnBasis_(input.decorrelationNumberOfIterations()),
    //lastValueB1_(0), estimatorB1_(0),nbValuesB1_(0),
    //estimatorB2_(input.decorrelationNumberOfIterations()),
    //lastValueD_(0), estimatorD_(0),nbValuesD_(0),
    //lastValueObservable_(0), meanObservable_(0),nbValuesObservable_(0),
    //lastValueBetterObservable_(0), meanBetterObservable_(0),nbValuesBetterObservable_(0),
    decorrelationNumberOfIterations_(input.decorrelationNumberOfIterations()),
    potential_(&potential)
  {}

  bool ControlVariate::isNone() const
  {
    return false;
  }
 
  double ControlVariate::potential(dvec const& position) const
  {
    return (*potential_)(position);
  }
  
  dvec ControlVariate::potentialDerivative(dvec const& position) const
  {
    return potential_->derivative(position);
  }
  
  double ControlVariate::potentialLaplacian(dvec const& position) const
  {
    return potential_->laplacian(position);
  }
  
  //### A MODIFIER
  void ControlVariate::appendToB(double observable, dvec const& position, dvec const& momentum, size_t indexOfIteration)
  {
    statsB_.append((observable - statsObservable_.mean()) * basisFunction(position, momentum), indexOfIteration);
    //lastValueB1_ = (observable - meanObservable_) * basisFunction(position, momentum);
    //estimatorB1_ = (nbValuesB1_* estimatorB1_ + lastValueB1_) / (nbValuesB1_+1);
    //nbValuesB1_++;
  }
  
  void ControlVariate::appendToD(double generatorOnBasisFunction, dvec const& position, dvec const& momentum, size_t indexOfIteration)
  {
    statsD_.append(- basisFunction(position, momentum) * generatorOnBasisFunction, indexOfIteration);
    
    //lastValueD_ = - basisFunction(position, momentum) * generatorOnBasisFunction;

    //estimatorD_ = (nbValuesD_* estimatorD_ + lastValueD_) / (nbValuesD_+1);
    //nbValuesD_++;
  }
  
  void ControlVariate::appendToObservable(double observable, size_t indexOfIteration)
  {
    statsObservable_.append(observable, indexOfIteration);
    
    /*lastValueObservable_ = observable;
    meanObservable_ = (nbValuesObservable_* meanObservable_ + lastValueObservable_) / (nbValuesObservable_+1);
    nbValuesObservable_++;*/
  }
  
  void ControlVariate::appendToBetterObservable(double observable, double generatorOnBasisFunction, dvec const& position, dvec const& momentum, size_t indexOfIteration)
  {
    if (statsD_.mean() != 0)
      statsBetterObservable_.append(observable + statsB_.mean() / statsD_.mean() * generatorOnBasisFunction, indexOfIteration);
    else
      statsBetterObservable_.append(observable, indexOfIteration);
    
    statsGeneratorOnBasis_.append(generatorOnBasisFunction, indexOfIteration);
    
    /*if (estimatorD_ != 0)
      lastValueBetterObservable_ = observable + statsB_.mean() / estimatorD_ * generatorOnBasisFunction;
    else
      lastValueBetterObservable_ = observable;
    
    //lastValueBetterObservable_ = generatorOnBasisFunction(position, momentum);
    
    meanBetterObservable_ = (nbValuesBetterObservable_* meanBetterObservable_ + lastValueBetterObservable_) / (nbValuesBetterObservable_+1);
    nbValuesBetterObservable_++;*/    
  }
  
  void ControlVariate::update(double observable, double generatorOnBasisFunction, dvec const& position, dvec const& momentum, size_t indexOfIteration)
  {
    //assert(generatorOnBasisFunction < 40000);
    appendToB(observable, position, momentum, indexOfIteration);
    appendToD(generatorOnBasisFunction, position, momentum, indexOfIteration);
    appendToObservable(observable, indexOfIteration);
    appendToBetterObservable(observable, generatorOnBasisFunction, position, momentum, indexOfIteration);
  }
  
  double ControlVariate::lastB() const
  {
    return statsB_.lastValue();
  }
  
  
  double ControlVariate::meanB() const
  {
    return statsB_.mean();
  }
  
  double ControlVariate::stdDeviationB() const
  {
    return statsB_.standardDeviation();
  }
  
  double ControlVariate::lastD() const
  {
    return statsD_.lastValue();
  }
  
  
  double ControlVariate::meanD() const
  {
    return statsD_.mean();
  }
  
  double ControlVariate::stdDeviationD() const
  {
    return statsD_.standardDeviation();
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
  
    double ControlVariate::lastGeneratorOnBasis() const
  {
    return statsGeneratorOnBasis_.lastValue();
  }
  
  
  double ControlVariate::meanGeneratorOnBasis() const
  {
    return statsGeneratorOnBasis_.mean();
  }
  
  double ControlVariate::stdDeviationGeneratorOnBasis() const
  {
    return statsGeneratorOnBasis_.standardDeviation();
  }
  
  double ControlVariate::autocorrelation(size_t indexDifference) const
  {
    return statsObservable_(indexDifference);
  }
  
  
  
  
  
  NoControlVariate::NoControlVariate(Input const& input, Potential& potential):
    ControlVariate(input, potential){}
    
  bool NoControlVariate::isNone() const
  {
    return true;
  }
    
  double NoControlVariate::basisFunction(dvec const& position, dvec const& momentum) const
  {
    cout << "No basisFunction for NoControlVariate" << endl;
    assert(false);
    return 0;
  }
  
    //double generatorOnBasisFunction(dvec const& position, dvec const& momentum) const;
  double NoControlVariate::laplacienQ(dvec const& position, dvec const& momentum) const
  {
    return 0;
  }
  
  dvec NoControlVariate::gradientQ(dvec const& position, dvec const& momentum) const
  {
    return dvec(position.size());
  }
  
  double NoControlVariate::laplacienP(dvec const& position, dvec const& momentum) const
  {
    return 0;
  }
    
  dvec NoControlVariate::gradientP(dvec const& position, dvec const& momentum) const
  {
    return dvec(position.size());
  }
  
  void NoControlVariate::update(double observable, double generatorOnBasisFunction, dvec const& position, dvec const& momentum, size_t indexOfIteration)
  {
    appendToObservable(observable, indexOfIteration);
  }

    
    
    
  
  
  
  SinusControlVariate::SinusControlVariate(Input const& input, Potential& potential):
    ControlVariate(input, potential)
  {}
  
  double SinusControlVariate::basisFunction(dvec const& position, dvec const& momentum) const
  {
    return sin(2 * M_PI * position(0));;
  }
   
  dvec SinusControlVariate::gradientQ(dvec const& position, dvec const& momentum) const
  {
    double q = position(0);
    dvec v(position.size());
    v(0) = 2 * M_PI * cos(2 * M_PI * position(0));
    return v;
  }
  
  double SinusControlVariate::laplacienQ(dvec const& position, dvec const& momentum) const
  {
    double q = position(0);
    return - pow (2 * M_PI, 2) * sin(2 * M_PI * position(0));
  }

  
  double SinusControlVariate::laplacienP(dvec const& position, dvec const& momentum) const
  {
    return 0;
  }
  
  dvec SinusControlVariate::gradientP(dvec const& position, dvec const& momentum) const
  {
    return dvec(position.size());
  }
  
  
  
  
    CosControlVariate::CosControlVariate(Input const& input, Potential& potential):
    ControlVariate(input, potential)
  {}
  
  double CosControlVariate::basisFunction(dvec const& position, dvec const& momentum) const
  {
    return cos(2 * M_PI * position(0));;
  }
  
  double CosControlVariate::laplacienQ(dvec const& position, dvec const& momentum) const
  {
    double q = position(0);
    return - pow (2 * M_PI, 2) * cos(2 * M_PI * position(0));
  }
  
  dvec CosControlVariate::gradientQ(dvec const& position, dvec const& momentum) const
  {
    double q = position(0);
    dvec v(position.size());
    v(0) = - 2 * M_PI * sin(2 * M_PI * position(0));
    return v;
  }
  
  double CosControlVariate::laplacienP(dvec const& position, dvec const& momentum) const
  {
    return 0;
  }
  
  dvec CosControlVariate::gradientP(dvec const& position, dvec const& momentum) const
  {
    return dvec(position.size());
  }
  
  
  
  
  
  
  
  
  
  SinExpControlVariate::SinExpControlVariate(Input const& input, Potential& potential):
    ControlVariate(input, potential)
  {}
  
  double SinExpControlVariate::basisFunction(dvec const& position, dvec const& momentum) const
  {
    return sin(2 * M_PI * position(0)) * exp(potential(position)/2);
  }
  
  dvec SinExpControlVariate::gradientQ(dvec const& position, dvec const& momentum) const
  {
    double q = position(0);
    dvec v(position.size());
    //v(0) = 2 * M_PI * cos(2 * M_PI * position(0));
    //v(0) = - 2 * M_PI * sin(2 * M_PI * position(0));
    v(0) = (2 * M_PI * cos(2 * M_PI * q)
	+ potentialDerivative(position)(0)/2 * sin(2 * M_PI * q))
	* exp(potential(position)/2);
    return v;
  }
  
    double SinExpControlVariate::laplacienQ(dvec const& position, dvec const& momentum) const
  {
    double q = position(0);
    
    return (-pow(2 * M_PI, 2) * sin(2* M_PI * q)
      +2 * M_PI * potentialDerivative(position)(0) *  cos(2 * M_PI * q)
      + potentialLaplacian(position)/2 * sin(2 * M_PI * q)
      + pow(potentialDerivative(position)(0) / 2, 2) * sin(2 * M_PI * q))
      * exp(potential(position)/2);
  }
  
  dvec SinExpControlVariate::gradientP(dvec const& position, dvec const& momentum) const
  {
    return dvec(position.size());
  }
    
  double SinExpControlVariate::laplacienP(dvec const& position, dvec const& momentum) const
  {
    return 0;
  }
  
  
  
  
  CosExpControlVariate::CosExpControlVariate(Input const& input, Potential& potential):
    ControlVariate(input, potential)
  {}
  
  double CosExpControlVariate::basisFunction(dvec const& position, dvec const& momentum) const
  {
    return cos(2 * M_PI * position(0)) * exp(potential(position)/2);
  }
  
  dvec CosExpControlVariate::gradientQ(dvec const& position, dvec const& momentum) const
  {
    double q = position(0);
    dvec v(position.size());
    //v(0) = 2 * M_PI * cos(2 * M_PI * position(0));
    //v(0) = - 2 * M_PI * sin(2 * M_PI * position(0));
    v(0) = (- 2 * M_PI * sin(2 * M_PI * q)
	+ potentialDerivative(position)(0)/2 * cos(2 * M_PI * q))
	* exp(potential(position)/2);
    return v;
  }
  
    double CosExpControlVariate::laplacienQ(dvec const& position, dvec const& momentum) const
  {
    double q = position(0);
    
    return (-pow(2 * M_PI, 2) * cos(2* M_PI * q)
      -2 * M_PI * potentialDerivative(position)(0) *  sin(2 * M_PI * q)
      + potentialLaplacian(position)/2 * cos(2 * M_PI * q)
      + pow(potentialDerivative(position)(0) / 2, 2) * cos(2 * M_PI * q))
      * exp(potential(position)/2);
  }
  
  dvec CosExpControlVariate::gradientP(dvec const& position, dvec const& momentum) const
  {
    return dvec(position.size());
  }
    
  double CosExpControlVariate::laplacienP(dvec const& position, dvec const& momentum) const
  {
    return 0;
  }
  
  
}