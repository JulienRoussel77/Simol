#include "simol/statphys/output/Statistics.hpp"

using std::min;
using std::max;

namespace simol
{
  
  //###### Statistics ######
  
  Statistics::Statistics(int nbRows, int nbCols):
    sumValues_(nbRows, nbCols),
    lastValue_(nbRows, nbCols),
    nbValues_(nbRows, nbCols),
    iidVar_(nbRows, nbCols)
  {
    sumValues_.fill(0);
    lastValue_.fill(0);
    nbValues_.fill(0);
    iidVar_.fill(0);
  };

  void Statistics::append(double value, int i, int j)
  {
    lastValue_(i, j) = value;
    sumValues_(i, j) += value;
    nbValues_(i, j)++;
  }


  double Statistics::mean(int i, int j) const
  {
    if (nbValues_(i, j) == 0)
      return 0;
    else
      return sumValues_(i, j) / nbValues(i, j);
  }

  Vector<double> Statistics::meanVec(int i) const
  {
    return piecewiseDivision(sumValues_.column(i), nbValues_.column(i));
  }

  DenseMatrix<double> Statistics::meanMat() const
  {
    return piecewiseDivision(sumValues_, nbValues_);
  }


  const long int& Statistics::nbValues(int i, int j) const
  {
    return nbValues_(i, j);
  }

  const Vector<long int> Statistics::nbValuesVec(int i) const
  {
    return nbValues_.column(i);
  }

  const DenseMatrix<long int>& Statistics::nbValuesMat() const
  {
    return nbValues_;
  }


  const double& Statistics::lastValue(int i, int j) const
  {
    return lastValue_(i, j);
  }

  const Vector<double> Statistics::lastValueVec(int i) const
  {
    return lastValue_.column(i);
  }

  const DenseMatrix<double>& Statistics::lastValueMat() const
  {
    return lastValue_;
  }  
  
  
  //###### IIDStats ######
  
  IIDStats::IIDStats(int nbRows, int nbCols):
    Statistics(nbRows, nbCols),
    iidVar_(nbRows, nbCols)
  {
    iidVar_.fill(0);
  };
  
  void IIDStats::append(double value, int i, int j)
  {
    lastValue_(i, j) = value;
    double tempTerm = iidVar_(i,j) + pow(mean(), 2);
    nbValues_(i, j)++;
    sumValues_(i, j) += value;
    iidVar_(i,j) = tempTerm - pow(mean(), 2) + (pow(value, 2) - tempTerm) / nbValues_(i, j);
  }
  
  const double& IIDStats::variance(int i, int j) const
  {
    return iidVar_(i, j);
  }
  
  double IIDStats::stdDeviation(int i, int j) const
  {
    return sqrt(variance(i,j));
  }
  
  //###### CorrelationStats ######


  CorrelationStats::CorrelationStats():
    decorrelationNbOfSteps_(0),
    timeStep_(0),
    nbOfAutocoPts_(0),
    nbOfObservables_(0),
    //timeStep_(timeStep),
    statsValues_(0),
    statsRefValues_(0),
    statsIntegratedCorrelation_(0),
    //statsIntegratedCorrelationTemp_(),
    statsCorrelation_(0, 0),
    indexRef_(0)
  {
    //cout << "CorrelationStats()" << endl;
  }


  CorrelationStats::CorrelationStats(int decorrelationNbOfSteps0, double timeStep0, int nbOfAutocoPts0, int nbOfObservables):
    decorrelationNbOfSteps_(decorrelationNbOfSteps0),
    timeStep_(timeStep0),
    nbOfAutocoPts_(nbOfAutocoPts0),
    nbOfObservables_(nbOfObservables),
    //timeStep_(timeStep),
    statsValues_(nbOfObservables),
    statsRefValues_(nbOfObservables),
    statsIntegratedCorrelation_(nbOfObservables),
    //statsIntegratedCorrelationTemp_(),
    statsCorrelation_(nbOfAutocoPts_, nbOfObservables),
    indexRef_(0),
    currentCorrelationIntegral_(0)
  {
    //cout << "CorrelationStats(int decorrelation : nbObs = " << decorrelationTime << endl;
  }
  
  const int& CorrelationStats::decorrelationNbOfSteps() const
  {
    return decorrelationNbOfSteps_;
  }
  
  double CorrelationStats::decorrelationTime() const
  {
    return decorrelationNbOfSteps() * timeStep();
  }
  
  const double& CorrelationStats::timeStep() const
  {
    return timeStep_;
  }
  
  void CorrelationStats::append(const double& newValue, long int iOfStep, int iOfObservable, const double& newRefValue)
  {
    //-- update the average --
    statsValues_.append(newValue, iOfObservable);
    //-- update the autocorrelations --
    if (decorrelationNbOfSteps_ != 0)
    {
      if (iOfStep % decorrelationNbOfSteps_ == 0)
      {
        statsRefValues_.append(newRefValue, iOfObservable);
        indexRef_ = iOfStep;
      }
      //-- update the correlation function --
      // decorrelationNbOfSteps_ / nbOfAutocoPts_ is the time intervale between two correlation estimations: can be bigger than the time step
      statsCorrelation_.append(statsRefValues_.lastValue(iOfObservable) * newValue, ((iOfStep - indexRef_)*nbOfAutocoPts_) / decorrelationNbOfSteps_, iOfObservable);
      //-- integration of the correlation function with a trapezoidal rule --
      currentCorrelationIntegral_ += ((indexRef_ == iOfStep) ? .5 : 1) * statsRefValues_.lastValue(iOfObservable) * newValue;
      if ((iOfStep+1) % decorrelationNbOfSteps_ == 0)    
      {
        statsIntegratedCorrelation_.append(timeStep() * currentCorrelationIntegral_, iOfObservable);
        currentCorrelationIntegral_ = 0;
      }
    }
  }

  double CorrelationStats::operator()(long int iOfStep, int iOfObservable) const
  {
    return statsCorrelation_.mean(iOfStep, iOfObservable);
  }


  const double& CorrelationStats::lastValue(int iOfObservable) const
  {
    return statsValues_.lastValue(iOfObservable);
  }


  double CorrelationStats::mean(int iOfObservable) const
  {
    return statsValues_.mean(iOfObservable);
  }


  int CorrelationStats::statsIntegratedCorrelationNbValues(int iOfObservable) const
  {
    return statsIntegratedCorrelation_.nbValues(iOfObservable);
  }


  double CorrelationStats::integratedCorrelation(int iOfObservable) const
  {
    return statsIntegratedCorrelation_.mean(iOfObservable);
  }

  Vector<double> CorrelationStats::integratedCorrelationVec() const
  {
    return statsIntegratedCorrelation_.meanMat().column(0);
  }

  double CorrelationStats::integratedCorrelationUnbiased(int iOfObservable) const
  {
    return integratedCorrelation(iOfObservable) - mean(iOfObservable) * statsRefValues_.mean(iOfObservable) * decorrelationTime();
  }
  
  double CorrelationStats::varIntegratedCorrelation(int iOfObservable) const
  {
    return statsIntegratedCorrelation_.variance(iOfObservable);
  }
  
  double CorrelationStats::varCorrelation(int indexDifference, int iOfObservable) const
  {
    return statsCorrelation_.variance(indexDifference, iOfObservable);
  }
  
  double CorrelationStats::stdDeviationCorrelation(int indexDifference, int iOfObservable) const
  {
    return sqrt(varCorrelation(indexDifference, iOfObservable));
  }
  
  double CorrelationStats::stdErrorCorrelation(int indexDifference, int iOfObservable) const
  {
    return stdDeviationCorrelation(indexDifference, iOfObservable) / sqrt(statsCorrelation_.nbValues(iOfObservable));
  }
  
  //###### AutocorrelationStats ######
  
  AutocorrelationStats::AutocorrelationStats():
    CorrelationStats()
  {}


  AutocorrelationStats::AutocorrelationStats(int decorrelationNbOfSteps0, double timeStep0, int nbOfAutocoPts0, int nbOfObservables):
    CorrelationStats(decorrelationNbOfSteps0, timeStep0, nbOfAutocoPts0, nbOfObservables)
  {}
  
  void AutocorrelationStats::append(double const& newValue, long int iOfStep, int iOfObservable)
  {
    CorrelationStats::append(newValue, iOfStep, iOfObservable, newValue);
  }

  double AutocorrelationStats::variance(int iOfObservable) const
  {
    return 2 * max(0., integratedCorrelationUnbiased(iOfObservable));
  }

  double AutocorrelationStats::standardDeviation(int iOfObservable) const
  {
    return sqrt(variance(iOfObservable));
  }
  
  ///
  ///Returns the variance of the estimator of the variance
  double AutocorrelationStats::varianceOfVariance(int iOfObservable) const
  {
    return 4 * varIntegratedCorrelation(iOfObservable);
  }
  
  ///
  ///Returns the variance of the estimator of the variance
  double AutocorrelationStats::stdDevOfVariance(int iOfObservable) const
  {
    return 2 * sqrt(varIntegratedCorrelation(iOfObservable));
  }
  
      ///
  ///Returns the variance of the estimator of the variance
  double AutocorrelationStats::stdErrorOfVariance(int iOfObservable) const
  {
    return stdDevOfVariance(iOfObservable) / sqrt(statsCorrelation_.nbValues(iOfObservable));
  }
  
}
