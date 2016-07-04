#include "simol/statphys/output/Statistics.hpp"

using std::min;
using std::max;

namespace simol
{
  Statistics::Statistics(int nbRows, int nbCols):
    sumValues_(nbRows, nbCols),
    lastValue_(nbRows, nbCols),
    nbValues_(nbRows, nbCols)
  {
    sumValues_.fill(0);
    lastValue_.fill(0);
    nbValues_.fill(0);
  };

  void Statistics::append(double value, int i, int j)
  {
    lastValue_(i, j) = value;
    if (nbValues_(i, j) == 0)
      sumValues_(i, j) = value;
    else
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


  const int& Statistics::nbValues(int i, int j) const
  {
    return nbValues_(i, j);
  }

  const Vector<int> Statistics::nbValuesVec(int i) const
  {
    return nbValues_.column(i);
  }

  const DenseMatrix<int>& Statistics::nbValuesMat() const
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


  AutocorrelationStats::AutocorrelationStats():
    decorrelationNbOfSteps_(0),
    decorrelationTime_(0),
    nbOfAutocoPts_(0),
    nbOfObservables_(0),
    //timeStep_(timeStep),
    statisticsValues_(0),
    statisticsRefValues_(0),
    statisticsMeanCorrelation_(0),
    //statisticsMeanCorrelationTemp_(),
    statisticsCorrelation_(0, 0),
    indexRef_(0)
  {
    //cout << "AutocorrelationStats()" << endl;
  }


  AutocorrelationStats::AutocorrelationStats(int decorrelationNbOfSteps0, double decorrelationTime, int nbOfAutocoPts0, int nbOfObservables):
    decorrelationNbOfSteps_(decorrelationNbOfSteps0),
    decorrelationTime_(decorrelationTime),
    nbOfAutocoPts_(nbOfAutocoPts0),
    nbOfObservables_(nbOfObservables),
    //timeStep_(timeStep),
    statisticsValues_(nbOfObservables),
    statisticsRefValues_(nbOfObservables),
    statisticsMeanCorrelation_(nbOfObservables),
    //statisticsMeanCorrelationTemp_(),
    statisticsCorrelation_(nbOfAutocoPts_, nbOfObservables),
    indexRef_(0)
  {
    //cout << "AutocorrelationStats(int decorrelation : nbObs = " << decorrelationTime << endl;
  }


  void AutocorrelationStats::append(double const& newValue, long int iOfStep, int iOfObservable)
  {
    append(newValue, iOfStep, iOfObservable, newValue);
  }


  double AutocorrelationStats::operator()(long int iOfStep, int iOfObservable) const
  {
    return statisticsCorrelation_.mean(iOfStep, iOfObservable);
  }


  const double& AutocorrelationStats::lastValue(int iOfObservable) const
  {
    return statisticsValues_.lastValue(iOfObservable);
  }


  double AutocorrelationStats::mean(int iOfObservable) const
  {
    return statisticsValues_.mean(iOfObservable);
  }


  int AutocorrelationStats::statisticsMeanCorrelationNbValues(int iOfObservable) const
  {
    return statisticsMeanCorrelation_.nbValues(iOfObservable);
  }


  double AutocorrelationStats::integratedAutocorrelation(int iOfObservable) const
  {
    return statisticsMeanCorrelation_.mean(iOfObservable) * decorrelationTime_;
  }


  Vector<double> AutocorrelationStats::integratedAutocorrelationVec() const
  {
    return statisticsMeanCorrelation_.meanMat().column(0) * decorrelationTime_;
  }




  double AutocorrelationStats::variance(int iOfObservable) const
  {
    return 2 * integratedAutocorrelationUnbiased(iOfObservable);
  }


  double AutocorrelationStats::standardDeviation(int iOfObservable) const
  {
    return sqrt(variance(iOfObservable));
  }

  void AutocorrelationStats::append(const double& newValue, long int iOfStep, int iOfObservable, const double& newRefValue)
  {
    //-- update the average --
    statisticsValues_.append(newValue, iOfObservable);
    //-- update the autocorrelations --
    if (decorrelationNbOfSteps_ != 0)
    {
      if (iOfStep % decorrelationNbOfSteps_ == 0)
      {
        statisticsRefValues_.append(newRefValue, iOfObservable);
        indexRef_ = iOfStep;
      }
      //-- update the correlation function --
      // decorrelationNbOfSteps_ / nbOfAutocoPts_ is the time intervale between two correlation estimations: can be bigger than the time step
      statisticsCorrelation_.append(statisticsRefValues_.lastValue(iOfObservable) * newValue, ((iOfStep - indexRef_)*nbOfAutocoPts_) / decorrelationNbOfSteps_, iOfObservable);
      //-- integration of the correlation function with a trapezoidal rule --
      statisticsMeanCorrelation_.append(((indexRef_ == iOfStep) ? .5 : 1) * statisticsRefValues_.lastValue(iOfObservable) * newValue, iOfObservable);
    }
  }

  double AutocorrelationStats::integratedAutocorrelationUnbiased(int iOfObservable) const
  {
    return max(0., integratedAutocorrelation(iOfObservable) - mean(iOfObservable) * statisticsRefValues_.mean(iOfObservable) * decorrelationTime_);
  }

}
