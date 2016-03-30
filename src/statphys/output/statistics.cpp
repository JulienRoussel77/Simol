#include "statistics.hpp"

using std::min;
using std::max;

namespace simol
{
  Statistics::Statistics(size_t nbRows, size_t nbCols):sumValues_(nbRows, nbCols), lastValue_(nbRows, nbCols), nbValues_(nbRows, nbCols)
  {
    //cout << "Statistics " << nbRows << " X " << nbCols << endl;
    sumValues_.fill(0);
    lastValue_.fill(0);
    nbValues_.fill(0);
  };

  void Statistics::append(double value, size_t i, size_t j)
  {
    lastValue_(i,j) = value;
    if (nbValues_(i,j) == 0)
      sumValues_(i,j) = value;
    else
      sumValues_(i,j) += value;
    nbValues_(i,j)++;
  }


  double Statistics::mean(size_t i, size_t j) const
  {
    if (nbValues_(i,j) == 0)
      return 0;
    else 
      return sumValues_(i,j) / nbValues(i,j);
  }
  
  Vector<double> Statistics::meanVec(size_t i) const
  {
    //return (sumValues_.array() / nbValues_.array().cast<double>()).matrix();
    return piecewiseDivision(sumValues_.column(i), nbValues_.column(i));
  }
  
  DenseMatrix<double> Statistics::meanMat() const
  {
    //return (sumValues_.array() / nbValues_.array().cast<double>()).matrix();
    return piecewiseDivision(sumValues_, nbValues_);
  }


  const size_t& Statistics::nbValues(size_t i, size_t j) const
  {
    return nbValues_(i,j);
  }

  const Vector<size_t> Statistics::nbValuesVec(size_t i) const
  {
    return nbValues_.column(i);
  }

  const DenseMatrix<size_t>& Statistics::nbValuesMat() const
  {
    return nbValues_;
  }


  const double& Statistics::lastValue(size_t i, size_t j) const
  {
    return lastValue_(i,j);
  }

  const Vector<double> Statistics::lastValueVec(size_t i) const
  {
    return lastValue_.column(i);
  }

  const DenseMatrix<double>& Statistics::lastValueMat() const
  {
    return lastValue_;
  }



  
  
  
  AutocorrelationStats::AutocorrelationStats():
    decorrelationNbOfIterations_(0),
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
  {//cout << "AutocorrelationStats()" << endl;
  }


  AutocorrelationStats::AutocorrelationStats(size_t decorrelationNbOfIterations0, double decorrelationTime, int nbOfAutocoPts0, size_t nbOfObservables):
    decorrelationNbOfIterations_(decorrelationNbOfIterations0),
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
  {//cout << "AutocorrelationStats(size_t decorrelation : nbObs = " << decorrelationTime << endl;
  }


  void AutocorrelationStats::append(double const& newValue, size_t iOfIteration, size_t iOfObservable)
  {
    append(newValue, iOfIteration, iOfObservable, newValue);
  }


  double AutocorrelationStats::operator()(size_t iOfIteration, size_t iOfObservable) const
  {
    return statisticsCorrelation_.mean(iOfIteration, iOfObservable);
  }


  const double& AutocorrelationStats::lastValue(size_t iOfObservable) const
  {
    return statisticsValues_.lastValue(iOfObservable);
  }


  double AutocorrelationStats::mean(size_t iOfObservable) const
  {
    return statisticsValues_.mean(iOfObservable);
  }


  size_t AutocorrelationStats::statisticsMeanCorrelationNbValues(size_t iOfObservable) const
  {
    return statisticsMeanCorrelation_.nbValues(iOfObservable);
  }


  double AutocorrelationStats::integratedAutocorrelation(size_t iOfObservable) const
  {
    return statisticsMeanCorrelation_.mean(iOfObservable) * decorrelationTime_;
  }


  Vector<double> AutocorrelationStats::integratedAutocorrelationVec() const
  {
    return statisticsMeanCorrelation_.meanMat().column(0) * decorrelationTime_;
  }




  double AutocorrelationStats::variance(size_t iOfObservable) const
  {
    return 2 * integratedAutocorrelationUnbiased(iOfObservable);
  }

  
  double AutocorrelationStats::standardDeviation(size_t iOfObservable) const
  {
    return sqrt(variance(iOfObservable));
  }
  
  void AutocorrelationStats::append(const double& newValue, size_t iOfIteration, size_t iOfObservable, const double& newRefValue)
  {
    if (iOfIteration % decorrelationNbOfIterations_ == 0)
    {
      statisticsRefValues_.append(newRefValue, iOfObservable);
      indexRef_ = iOfIteration;
    }
    statisticsValues_.append(newValue, iOfObservable);
		//On intègre le profil d'autocorrélation à l'ordre 2 en coomptant la valeur en 0 pour un demi
    statisticsMeanCorrelation_.append(((indexRef_ == iOfIteration)?.5:1) * statisticsRefValues_.lastValue(iOfObservable) * newValue, iOfObservable);
    statisticsCorrelation_.append(statisticsRefValues_.lastValue(iOfObservable) * newValue, ((iOfIteration - indexRef_)*nbOfAutocoPts_) / decorrelationNbOfIterations_, iOfObservable);    
  }
  
  /*template <>
  void AutocorrelationStats<Vector<double>>::append(Vector<double> const& newValue, size_t iOfIteration, size_t iOfObservable, Vector<double> const& newRefValue)
  {
    if (iOfIteration % decorrelationNbOfIterations_ == 0)
    {
      statisticsRefValues_.append(newValue, iOfObservable);
      indexRef_ = iOfIteration;
    }
    statisticsValues_.append(newValue, iOfObservable);
		//On intègre le profil d'autocorrélation à l'ordre 2
    statisticsMeanCorrelation_.append(((indexRef_ == iOfIteration)?.5:1) * dot(statisticsRefValues_.lastValue(iOfObservable), newValue), iOfObservable);
    statisticsCorrelation_.append(dot(statisticsRefValues_.lastValue(iOfObservable), newValue), ((iOfIteration - indexRef_)*nbOfAutocoPts_) / decorrelationNbOfIterations_, iOfObservable);
  }*/
  
  double AutocorrelationStats::integratedAutocorrelationUnbiased(size_t iOfObservable) const
  {
    return max(0., integratedAutocorrelation(iOfObservable) - mean(iOfObservable) * statisticsRefValues_.mean(iOfObservable) * decorrelationTime_);
  }
  
  /*double AutocorrelationStats<Vector<double>>::integratedAutocorrelationUnbiased(size_t iOfObservable) const
  {
    return max(0., integratedAutocorrelation(iOfObservable) - dot(mean(iOfObservable), statisticsRefValues_.mean(iOfObservable)) * decorrelationTime_);
  }*/
  
}
