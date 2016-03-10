#include "statistics.hpp"

namespace simol{

  template <>
  void AutocorrelationStats<double>::append(double const& newValue, size_t iOfIteration, size_t iOfObservable, double const& newRefValue)
  {
    if (iOfIteration % decorrelationNumberOfIterations_ == 0)
    {
      statisticsRefValues_.append(newRefValue, iOfObservable);
      indexRef_ = iOfIteration;
      //cout << iOfObservable << " : " << statisticsRefValues_.lastValue(iOfObservable) << " * " << newValue << " = " << statisticsRefValues_.lastValue(iOfObservable) * newValue << " __" <<newRefValue << endl;

    }
    statisticsValues_.append(newValue, iOfObservable);
    //cout << iOfObservable << " : " << statisticsRefValues_.lastValue(iOfObservable) << " * " << newValue << " = " << statisticsRefValues_.lastValue(iOfObservable) * newValue << endl;
    statisticsMeanCorrelation_.append(statisticsRefValues_.lastValue(iOfObservable) * newValue, iOfObservable);
    statisticsCorrelation_.append(statisticsRefValues_.lastValue(iOfObservable) * newValue, iOfIteration - indexRef_, iOfObservable);
  }

  template <>
  void AutocorrelationStats<Vector<double>>::append(Vector<double> const& newValue, size_t iOfIteration, size_t iOfObservable, Vector<double> const& /*newRefValue*/)
  {
    if (iOfIteration % decorrelationNumberOfIterations_ == 0)
    {
      statisticsRefValues_.append(newValue, iOfObservable);
      indexRef_ = iOfIteration;
      /*if (iOfIteration > 0)
      {
	statisticsMeanCorrelation_.append(statisticsMeanCorrelationTemp_.mean());
	statisticsMeanCorrelationTemp_.clear();
      }*/
    }
    statisticsValues_.append(newValue, iOfObservable);
    //statisticsMeanCorrelationTemp_.append(dot(statisticsRefValues_.lastValue(), newValue));
    statisticsMeanCorrelation_.append(dot(statisticsRefValues_.lastValue(iOfObservable), newValue), iOfObservable);
    statisticsCorrelation_.append(dot(statisticsRefValues_.lastValue(iOfObservable), newValue), iOfIteration - indexRef_, iOfObservable);
  }

  template <>
  double AutocorrelationStats<double>::integratedAutocorrelationUnbiased(size_t iOfObservable) const
  {
    return std::max(0., integratedAutocorrelation(iOfObservable) - mean(iOfObservable) * statisticsRefValues_.mean(iOfObservable) * decorrelationTime_);
  }

  template <>
  double AutocorrelationStats<Vector<double>>::integratedAutocorrelationUnbiased(size_t iOfObservable) const
  {
    return std::max(0., integratedAutocorrelation(iOfObservable) - dot(mean(iOfObservable), statisticsRefValues_.mean(iOfObservable)) * decorrelationTime_);
  }

}
