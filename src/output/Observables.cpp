#include "simol/output/Observable.hpp"

using std::cout;
using std::endl;
using std::sin;
using std::cos;
using std::exp;
using std::ofstream;

namespace simol
{
  
  /*Observable* createObservable(const Input& input)
  {
    return new Observable(input);
  }*/

  Observable::Observable(Input const& input, int idObs, int decorrelationNbOfSteps0, int nbOfAutocoPts0):
    decorrelationNbOfSteps_(decorrelationNbOfSteps0),
    timeStep_(input.timeStep()),
    printPeriodNbOfSteps_(input.printPeriodNbOfSteps()),
    nbOfAutocoPts_(nbOfAutocoPts0),
    currentValue_(0),
    autocoStats_(decorrelationNbOfSteps(), timeStep(), nbOfAutocoPts()),
    outPath_(input.nameOfObs(idObs)+".txt"),
    outFlux_(input.outputFolderName() + input.nameOfObs(idObs)+".txt")
  {
    outFlux() << "# time value mean variance" << endl;
    //if (outPathFinal != "")
    //  outFluxFinal_ = std::make_shared<ofstream>(input.simuTypeName() + outPathFinal, std::ofstream::app);
    //cout << "-----   " << input.simuTypeName() + outPathFinal << endl;
    
    if (doComputeCorrelations())
    {
      outFluxCorrelation_ = make_shared<ofstream>(input.outputFolderName() + "correlation_" + input.nameOfObs(idObs)+".txt");
      outFluxCorrelation() << "# iOfSpan correlation varOfCorrelation" << endl;
    }
  }
  
  Observable::~Observable(){}
  
  

  const double& Observable::timeStep() const
  {
    return timeStep_;
  }

  const int& Observable::decorrelationNbOfSteps() const
  {
    return decorrelationNbOfSteps_;
  }

  double Observable::decorrelationTime() const
  {
    return decorrelationNbOfSteps() * timeStep();
  }

  int const& Observable::nbOfAutocoPts() const
  {return nbOfAutocoPts_;}
  
  double Observable::autocoPtsPeriod() const
  {return decorrelationTime() / nbOfAutocoPts();}
  
  double Observable::printPeriodTime() const
  {return printPeriodNbOfSteps_ * timeStep();}

  const int& Observable::printPeriodNbOfSteps() const
  {return printPeriodNbOfSteps_;}
  
  double Observable::time() const
  {return autocoStats_.nbValues() * printPeriodTime();}
  
  string const& Observable::outPath() const
  {return outPath_;}
  
  ofstream& Observable::outFlux()
  {
    return outFlux_;
  }
  
  ofstream& Observable::outFluxCorrelation()
  {
    return *outFluxCorrelation_;
  }
  
  bool Observable::doComputeCorrelations() const
  {return decorrelationNbOfSteps()>0;}


  const double& Observable::currentValue() const
  {return currentValue_;}
  
  double& Observable::currentValue()
  {return currentValue_;}

  double Observable::lastValue() const
  {
    return autocoStats_.lastValue();
  }


  double Observable::mean() const
  {
    return autocoStats_.mean();
  }
  
  double Observable::asymptoticVariance() const
  {
    return autocoStats_.asymptoticVariance();
  }
  
  double Observable::asyvarOfAsyvar() const
  {
    return autocoStats_.asyvarOfAsyvar();
  }  

  double Observable::correlationAtSpan(int iOfSpan) const
  {
    return autocoStats_.correlationAtSpan(iOfSpan);
  }
  
  double Observable::centeredCorrelationAtSpan(int iOfSpan) const
  {
    //return autocoStats_(iOfSpan) - pow(meanObservable(), 2);
    return autocoStats_.centeredCorrelationAtSpan(iOfSpan);
  }
  
  double Observable::varCorrelationAtSpan(int iOfSpan) const
  {
    return autocoStats_.varCorrelationAtSpan(iOfSpan);
  }

  double Observable::stdDevCorrelationAtSpan(int iOfSpan) const
  {
    return autocoStats_.stdDevCorrelationAtSpan(iOfSpan);
  }
  
  
  
  void Observable::append(double value, long int iOfStep)
  {
    autocoStats_.append(value, iOfStep);
  }
  
  void Observable::appendCurrent(long int iOfStep)
  {
    append(currentValue(), iOfStep);
  }
  
  
  bool Observable::doOutput(long int iOfStep) const
  {
    return (printPeriodNbOfSteps_ > 0 && iOfStep % printPeriodNbOfSteps_ == 0);
  }
  
  /*bool Observable::doFinal() const
  {
    return outFluxFinal_ != nullptr;
  }*/


  void Observable::display(long int iOfStep)
  {
    //cout << "10play !" << endl;
    outFlux() << iOfStep * timeStep() 
             << " " << lastValue()
             << " " << mean()
             << " " << asymptoticVariance()
             << " " << asyvarOfAsyvar() << endl;
             //<< " " << autocoStats_.integratedCorrelation()<< endl;
  }
  
  void Observable::displayFinalValues(ofstream& out)
  {
    out << mean()
        << " " << asymptoticVariance()
        << " " << asyvarOfAsyvar() << endl;
  }
  
  void Observable::displayCorrelations(long int iOfStep)
  {
    //cout << "Velocity : The correlation in 0 is " << floor((2 * obsVelocity().centeredCorrelationAtSpan(0) * decorrelationTime() / nbOfAutocoPts()) / obsVelocity().asymptoticVariance() * 10000)/100 << "% of the variance" << endl;
    //cout << "SumFlux : The correlation in 0 is " << floor((2 * obsSumFlux().centeredCorrelationAtSpan(0) * decorrelationTime() / nbOfAutocoPts()) / obsSumFlux().asymptoticVariance() * 10000)/100 << "% of the variance" << endl;
    //cout << "ModiFlux : The correlation in 0 is " << floor((2 * obsModiFlux().centeredCorrelationAtSpan(0) * decorrelationTime() / nbOfAutocoPts()) / obsModiFlux().asymptoticVariance() * 10000)/100 << "% of the variance" << endl;
       
    //cout << "Velocity : The correlation in 0 is " << floor((2 * obsVelocity().centeredCorrelationAtSpan(0) * decorrelationTime() / nbOfAutocoPts()) / obsVelocity().asymptoticVariance() * 10000)/100 << "% of the variance" << endl;
    outFluxCorrelation() << "# The correlation in iOfSpan = 0 is " << floor((2 * centeredCorrelationAtSpan(0) * decorrelationTime() / nbOfAutocoPts()) / asymptoticVariance() * 10000)/100 << "% of the variance" << endl;
    
    for (int iOfSpan = 0; iOfSpan < nbOfAutocoPts(); iOfSpan++)
    {
      //cout << sqrt(varCorrelationAtSpan(iOfSpan)) << " / " << sqrt(iOfStep * timeStep()) << " = " << sqrt(varCorrelationAtSpan(iOfSpan) / (iOfStep * timeStep())) << endl;
      outFluxCorrelation() << iOfSpan * autocoPtsPeriod()
                       << " " << centeredCorrelationAtSpan(iOfSpan)
                       << " " << sqrt(varCorrelationAtSpan(iOfSpan) / (iOfStep * timeStep())) << endl;
    }
  }
}










