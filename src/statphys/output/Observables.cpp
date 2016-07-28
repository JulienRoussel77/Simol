#include "simol/statphys/output/Observable.hpp"

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

  Observable::Observable(Input const& input, const string& outPath):
    decorrelationNbOfSteps_(input.decorrelationNbOfSteps()),
    timeStep_(input.timeStep()),
    printPeriodNbOfSteps_(input.printPeriodNbOfSteps()),
    nbOfAutocoPts_(input.nbOfAutocoPts()),
    currentValue_(0),
    autocoStats_(decorrelationNbOfSteps(), timeStep(), nbOfAutocoPts()),
    outFlux_(input.outputFolderName() + outPath)
  {
    outFlux() << "# time value mean variance" << endl;
    //if (outPathFinal != "")
    //  outFluxFinal_ = std::make_shared<ofstream>(input.simuTypeName() + outPathFinal, std::ofstream::app);
    //cout << "-----   " << input.simuTypeName() + outPathFinal << endl;
    
    if (doComputeCorrelations())
    {
      outFluxCorrelation_ = new ofstream(input.outputFolderName() + "correlation_" + outPath);
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
  
  double Observable::variance() const
  {
    return autocoStats_.variance();
  }

  double Observable::stdDev() const
  {
    return autocoStats_.stdDev();
  }
  
  double Observable::varOfVar() const
  {
    return autocoStats_.varOfVar();
  }
  
  double Observable::stdDevOfVar() const
  {
    return autocoStats_.stdDevOfVar();
  }
  
  

  double Observable::correlationAtSpan(int iOfSpan) const
  {
    return autocoStats_.correlationAtSpan(iOfSpan);
  }
  
  double Observable::unbiasedCorrelationAtSpan(int iOfSpan) const
  {
    //return autocoStats_(iOfSpan) - pow(meanObservable(), 2);
    return autocoStats_.unbiasedCorrelationAtSpan(iOfSpan);
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
    currentValue() = 0;
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
    outFlux() << iOfStep * timeStep() 
             << " " << lastValue()
             << " " << mean()
             << " " << variance()
             << " " << varOfVar() << endl;
  }
  
  /*void Observable::finalDisplay(vector<double> parameters)
  {
    if (doFinal())
    {
      for (int iOfParameter=0; iOfParameter < (int)parameters.size(); iOfParameter++)
        *outFluxFinal_  << parameters[iOfParameter] << " " << setw(6);
        *outFluxFinal_ << " " << setw(12) << mean()
                       << " " << setw(12) << variance()
                       << " " << setw(12) << varOfVar() << endl;
    }
  }*/
  
  void Observable::displayCorrelations(long int iOfStep)
  {
    //cout << "Velocity : The correlation in 0 is " << floor((2 * obsVelocity().unbiasedCorrelationAtSpan(0) * decorrelationTime() / nbOfAutocoPts()) / obsVelocity().variance() * 10000)/100 << "% of the variance" << endl;
    //cout << "SumFlow : The correlation in 0 is " << floor((2 * obsSumFlow().unbiasedCorrelationAtSpan(0) * decorrelationTime() / nbOfAutocoPts()) / obsSumFlow().variance() * 10000)/100 << "% of the variance" << endl;
    //cout << "ModiFlow : The correlation in 0 is " << floor((2 * obsModiFlow().unbiasedCorrelationAtSpan(0) * decorrelationTime() / nbOfAutocoPts()) / obsModiFlow().variance() * 10000)/100 << "% of the variance" << endl;
       
    //cout << "Velocity : The correlation in 0 is " << floor((2 * obsVelocity().unbiasedCorrelationAtSpan(0) * decorrelationTime() / nbOfAutocoPts()) / obsVelocity().variance() * 10000)/100 << "% of the variance" << endl;
    outFluxCorrelation() << "# The correlation in iOfSpan = 0 is " << floor((2 * unbiasedCorrelationAtSpan(0) * decorrelationTime() / nbOfAutocoPts()) / variance() * 10000)/100 << "% of the variance" << endl;
    
    for (int iOfSpan = 0; iOfSpan < nbOfAutocoPts(); iOfSpan++)
    {
      outFluxCorrelation() << iOfSpan * autocoPtsPeriod()
                       << " " << unbiasedCorrelationAtSpan(iOfSpan)
                       << " " << sqrt(varCorrelationAtSpan(iOfSpan) / (iOfStep * timeStep())) << endl;
    }
  }
}










