#include "statistics.hpp"

using std::min;
using std::max;


namespace simol{  
  
  template <>
  void AutocorrelationStats<double>::append(double const& newValue, int indexOfIteration)
  {
    if (indexOfIteration % decorrelationNumberOfIterations_ == 0)
    {
      statisticsRefValues_.append(newValue);
      indexRef_ = indexOfIteration;
    }
    statisticsCorrelation_.append(statisticsRefValues_.lastValue() * newValue, indexOfIteration - indexRef_);
    statisticsValues_.append(newValue);
    statisticsMeanCorrelation_.append(statisticsRefValues_.lastValue() * newValue);
  }
  
  template <>
  void AutocorrelationStats<dvec>::append(dvec const& newValue, int indexOfIteration)
  {
    if (indexOfIteration % decorrelationNumberOfIterations_ == 0)
    {
      statisticsRefValues_.append(newValue);
      indexRef_ = indexOfIteration;
      if (indexOfIteration > 0)
	statisticsMeanCorrelation_.append(statisticsMeanCorrelationTemp_.mean());
    }
    statisticsValues_.append(newValue);
    /*cout << "newValue : " << newValue 
    << " / nbVal = " << statisticsValues_.nbValues(0)
    << " / mean = " << statisticsValues_() << endl;*/
    statisticsMeanCorrelationTemp_.append(statisticsRefValues_.lastValue().dot(newValue));
    statisticsCorrelation_.append(statisticsRefValues_.lastValue().dot(newValue), indexOfIteration - indexRef_);
    //cout << indexOfIteration - indexRef_ << " / " << indexOfIteration << " -> " << valueRef_.dot(newValue) << endl;
  }
  
  template <>
  double AutocorrelationStats<double>::integratedAutocorrelationUnbiased() const
  {
    return max(0., integratedAutocorrelation() - mean() * statisticsRefValues_.mean());// * decorrelationNumberOfIterations_ * timeStep());
  }
  
  template <>
  double AutocorrelationStats<dvec>::integratedAutocorrelationUnbiased() const
  {
    return max(0., integratedAutocorrelation() - mean().dot(statisticsRefValues_.mean()));// * decorrelationNumberOfIterations_ * timeStep());
  }
  
}